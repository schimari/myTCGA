# analyses of TCGA data retrieved from r_retrieve.r

library("survival")
library("limma")
library("dplyr")
library("survminer")
library("RColorBrewer")

cohort <- "KIRC"

# read RNA file 
rna <- read.table(file="KIRC/2016-01-28/rnaseq/KIRC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=T,row.names=1,sep="\t")
# and take off first row cause we don’t need it
rna <- rna[-1,]
# and read the Clinical file, in this case I transposed it to keep the clinical feature title as column name
clinical <- t(read.table(file="KIRC/2016-01-28/clinical/KIRC.Merge_Clinical.Level_1/KIRC.merged_only_clinical_clin_format.txt", header=T, row.names=1, sep="\t"))

# color
colfunc <- colorRampPalette(rev(c("darkred", "darkorange", "darkgreen", "darkblue")))

# first I remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]

#now I need to identify normal and tumor samples. this is done using the TCGA barcode (https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode). 
#The two digits at position 14-15 of the barcode will indicate the sample type, from the link:
#"Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29."

# see the values
table(substr(colnames(rna),14,14))

# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == "1")	# normal
t_index <- which(substr(colnames(rna),14,14) == "0")	# tumor

# apply voom function from limma package to normalize the data
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}

rna_vm  <- vm(rna)
colnames(rna_vm) <- gsub("\\.","-",substr(colnames(rna),1,12))

# and check how data look, they should look normally-ish distributed
#hist(rna_vm)
# we can remove the old “rna” cause we don’t need it anymor
rm(rna)

# calculate z-scores
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])
# set the rownames keeping only gene name
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,"\\|"))[[1]])

rm(rna_vm) #we don’t need it anymore

# match the patient ID in clinical data with the colnames of z_rna
clinical <- as.data.frame(clinical)
clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
sum(clinical$IDs %in% colnames(z_rna)) # we have 533 patients that we could use

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
ind_keep <- grep("days_to_new_tumor_event_after_initial_treatment",colnames(clinical))

# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if(sum(is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- max(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,"NA")
  }
}

# do the same to death
ind_keep <- grep("days_to_death",colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if(sum(is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,"NA")
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep("days_to_last_followup",colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if(sum(is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,"NA")
  }
}

# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")

# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                 as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

# create vector for death censoring
table(clinical$patient.vital_status)
# alive  dead 
#   372   161 

all_clin$death_event <- ifelse(clinical$patient.vital_status == "alive", 0,1)

#finally add row.names to clinical
rownames(all_clin) <- clinical$IDs

# create event vector for RNASeq data
#event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))	# only up regulated
event_rna <- t(apply(z_rna, 1, function(x) ifelse(x > 1.96,1,ifelse(x < -1.96,2,0))))

# since we need the same number of patients in both clinical and RNASeq data take the indices for the matching samples
ind_tum <- which(unique(colnames(z_rna)) %in% rownames(all_clin))
ind_clin <- which(rownames(all_clin) %in% colnames(z_rna))

# pick your gene of interest
geneint <- "CDK11B"
ind_gene <- which(rownames(z_rna) == geneint)

# check how many altered samples we have
table(event_rna[ind_gene,])

# run survival analysis
s <- survfit(Surv(as.numeric(as.character((all_clin$new_death/365)*12))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum])
#s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

# extraect the p.value
#pv <- ifelse(is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),6)))[[1]]

# survminer plot
cairo_pdf(sprintf("plots/%s_%s.pdf", geneint, cohort), width=11.69, height=8.27)
ggsurvplot(s,  size = 0.6,  # change line size
	#main="",
	palette = colfunc(nrow(s)), # custom color palettes
	conf.int = TRUE, # Add confidence interval
	pval = TRUE, # Add p-value
	risk.table = TRUE, # Add risk table
	risk.table.col = "strata", # Risk table color by groups
	legend.labs = c(sprintf("%s not altered", geneint), sprintf("%s up", geneint), sprintf("%s down", geneint)), # Change legend labels
	risk.table.height = 0.25, # Useful to change when you have multiple groups
	surv.scale = "percent",
	xlab = "Time (months)",
	ggtheme = theme_bw()) # Change ggplot2 theme
dev.off()
	   





