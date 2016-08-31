## TCGA with limma and topTreat

library("survival")
library("limma")
library("RColorBrewer")
library("dplyr")

cohort <- "KIRC"
fdate <- "2016-01-28"
bfolder <- sprintf("%s/%s/d_vs_a", cohort, fdate)

# read RNA file 
rna <- read.table(file=sprintf("%s/%s/data/rnaseq/KIRC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", cohort, fdate), header=T,row.names=1,sep="\t")
# and take off first row cause we don’t need it
rna <- rna[-1,]
# and read the Clinical file, in this case I transposed it to keep the clinical feature title as column name
clinical <- t(read.table(file=sprintf("%s/%s/data/clinical/KIRC.Merge_Clinical.Level_1/KIRC.merged_only_clinical_clin_format.txt", cohort, fdate), header=T, row.names=1, sep="\t"))

# color
colfunc <- colorRampPalette(rev(c("darkred", "darkorange", "darkgreen", "darkblue")))

# folder stuff
if(!dir.exists(bfolder)) { 
	dir.create(bfolder)
}

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

# now I need to identify normal and tumor samples. this is done using the TCGA barcode (https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode). 
# The two digits at position 14-15 of the barcode will indicate the sample type, from the link:
# "Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29."

# see the values
table(substr(colnames(rna),14,14))

# get the index of the tumor/control samples
n_index <- which(substr(colnames(rna),14,14) == "1")	# normal
t_index <- which(substr(colnames(rna),14,14) == "0")	# tumor

# apply voom function from limma package to normalize the data
# Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship 
# and use this to compute appropriate observational-level weights. The data are then ready for linear modelling.
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F,normalize.method = "none")	# or plot=T for a mean-variance trend plot
  return(ex$E)
}

rna_vm  <- vm(rna)
colnames(rna_vm) <- gsub("\\.","-",substr(colnames(rna),1,12))
rownames(rna_vm) <- sapply(rownames(rna_vm), function(x) unlist(strsplit(x,"\\|"))[[1]])

# select only tumors
t_rna_vm <- rna_vm[,t_index]
t_col <- colnames(t_rna_vm)

# get corresponding clinical data
t_clin_all <- clinical[clinical$IDs %in% t_col,]

# set dead or alive
t_clin_alive <- t_clin_all[t_clin_all$patient.vital_status=="alive",]
t_clin_dead <- t_clin_all[t_clin_all$patient.vital_status=="dead",]

t_ca <- t_clin_alive$IDs
t_cd <- t_clin_dead$IDs

# get corresponding RNA counts
t_rna_alive <- as.data.frame(t_rna_vm[, colnames(t_rna_vm) %in% t_ca])
t_rna_dead <- as.data.frame(t_rna_vm[, colnames(t_rna_vm) %in% t_cd])

# make unique row names before cbind
nams <- rownames(t_rna_alive)
rownames(t_rna_alive) <- make.names(nams, unique=TRUE)
rownames(t_rna_dead) <- make.names(nams, unique=TRUE)

t_rna_df <- cbind(t_rna_alive, t_rna_dead)

# prepare model matrix
x <- NULL
for (i in 1:length(colnames(t_rna_alive))) { 
	x <- append(x, "1") 
	}
for (i in 1:length(colnames(t_rna_dead))) { 
	x <- append(x, "2") 
	}

Faktoren <- as.numeric(x)
design <- model.matrix(~ 0+factor(Faktoren))
colnames(design) <- c("alive", "dead")

# do linear fitting
fit <- lmFit(t_rna_df, design)
fit <- eBayes(fit)

# contrast between dead and alive
contrast.matrix <- makeContrasts(dead-alive, levels=design)
fitCon <- contrasts.fit(fit, contrast.matrix)
fitTreat <- treat(fitCon)
fitCon <- eBayes(fitCon)

# toptreat, etc.
t_treat <- topTreat(fitTreat, coef=1, number=20, adjust.method = "BH")
t_logFC <- toptable(fitCon, coef=1, number=20, sort.by="logFC", adjust.method = "BH")

# cbioPortal 
# logFC Genes :	http://www.cbioportal.org/index.do?cancer_study_list=kirc_tcga_pub&cancer_study_id=kirc_tcga_pub&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=kirc_tcga_pub_mutations&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=kirc_tcga_pub_gistic&genetic_profile_ids_PROFILE_MRNA_EXPRESSION=kirc_tcga_pub_rna_seq_v2_mrna_median_Zscores&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=kirc_tcga_pub_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list=SAA1%0D%0ASAA2+++%0D%0APPP1R1A+%0D%0APLG+++++%0D%0APAEP++++%0D%0ASLC6A19+%0D%0ATNNT1+++%0D%0AIL20RB++%0D%0ALBP+++++%0D%0AMAT1A+++%0D%0APTPRH+++%0D%0AHP++++++%0D%0APITX1+++%0D%0AIGFN1+++%0D%0ACOL7A1++%0D%0ASLPI++++%0D%0AF2++++++%0D%0AUGT1A10+%0D%0ASLC16A12%0D%0ASERPINA3%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit&show_samples=false&
# topTreat Genes : http://www.cbioportal.org/index.do?cancer_study_list=kirc_tcga_pub&cancer_study_id=kirc_tcga_pub&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=kirc_tcga_pub_mutations&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=kirc_tcga_pub_gistic&genetic_profile_ids_PROFILE_MRNA_EXPRESSION=kirc_tcga_pub_rna_seq_v2_mrna_median_Zscores&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=kirc_tcga_pub_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list=ZIC2+%0D%0ACDCA3+++%0D%0AADAMTS14%0D%0ASOWAHB+%0D%0AFKBP11++%0D%0ARHNO1%0D%0AHN1+++++%0D%0ATNNT1+++%0D%0AB3GNTL1+%0D%0AUBE2C+++%0D%0AWDR72+++%0D%0AAMOT++++%0D%0ASHOX2+++%0D%0ATROAP+++%0D%0ANOP2++++%0D%0ANUMBL+++%0D%0AITGA6+++%0D%0AITPKA+++%0D%0AWFDC10B+%0D%0APDCD5+++%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit&show_samples=false&











