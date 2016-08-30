## TCGA with limma and topTreat

library("survival")
library("limma")
library("RColorBrewer")
library("dplyr")

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
  ex <- voom(x,d,plot=F)	# or plot=T for a mean-variance trend plot
  #ex <- voom(x,d,plot=F,normalize.method = "none")	# or plot=T for a mean-variance trend plot
  return(ex$E)
}

rna_vm  <- vm(rna)
colnames(rna_vm) <- gsub("\\.","-",substr(colnames(rna),1,12))
rownames(rna_vm) <- sapply(rownames(rna_vm), function(x) unlist(strsplit(x,"\\|"))[[1]])

## linear modelling
# 1 = Tumor, 2 = Normal
tmpdf1 <- data.frame(t=t_index, matrix=1)
tmpdf2 <- data.frame(t=n_index, matrix=2)
tmpdf <- rbind(tmpdf1, tmpdf2)
tmpdfs <- tmpdf[sort.list(tmpdf$t),]

Faktoren <- as.numeric(tmpdfs$matrix)
design <- model.matrix(~ 0+factor(Faktoren))
colnames(design) <- c("tumor", "normal")

fit <- lmFit(rna_vm, design)
fit <- eBayes(fit)	# not sure if it is necessary

# contrasts
contrast.matrix <- makeContrasts(tumor-normal, levels=design)
fitCon <- contrasts.fit(fit, contrast.matrix)
fitTreat <- treat(fitCon, lfc=.5)	#or lfc=log2(1.2)
fitCon <- eBayes(fitCon)

ttlog <- toptable(fitCon, coef=1, number=20, sort.by="logFC", adjust.method = "BH")
colnames(ttlog) <- c("Symbols","logFC","t","P.Value","adj.P.Val","B")
ttpval <- toptable(fitCon, coef=1, number=20, adjust.method = "BH")
colnames(ttpval) <- c("Symbols","logFC","t","P.Value","adj.P.Val","B")
# I prefer topTreat
tttreat <- topTreat(fitTreat, coef=1, number=20, adjust.method = "BH")
colnames(tttreat) <- c("Symbols","logFC","AveExpr","t","P.Value","adj.P.Val")

# create genelist for survival plots with x_TCGA.r
genelist <- ""
for(i in 1:length(tttreat$Symbols)){
	genelist <- paste(genelist, tttreat$Symbols[i], sep="\n")
}
write.table(genelist, file="KIRC/x_tables/genelist_topTreat.txt", quote=F, col.names=F, row.names=F)

# plotting
filename <- "KIRC/x_plots/prad_top_genes.pdf"
cairo_pdf(filename, width=8.27, height=11.69)
par(mfrow=c(3,1), cex.lab=1.5, mar=c(7, 5, 4, 2) + 0.1)	# mar: c(bottom, left, top, right)
barplot(ttlog$logFC, beside=T, main="KIRC most differentially expressed Genes ordered by log2(fold-change)", names.arg=ttlog$Symbols, col=ifelse(ttlog$logFC>0,"darkgreen","darkred"), las=3, cex.names = 1.2, cex.main = 1.2, axis.lty = 1, axes=T, ylab="logFC", ylim = c(-10,10))
abline(h=0)
barplot(ttpval$logFC, beside=T, main="KIRC most differentially expressed Genes ordered by adjusted p-value (BH)", names.arg=ttpval$Symbols, col=ifelse(ttpval$logFC>0,"darkgreen","darkred"), las=3, cex.names = 1.2, cex.main = 1.2, axis.lty = 1, axes=T, ylab="logFC", ylim = c(-10,10))
abline(h=0)
barplot(tttreat$logFC, beside=T, main="KIRC most differentially expressed Genes (treat) ordered by adjusted p-value (BH)\nlog(fold-change) TH = 0.5", names.arg=tttreat$Symbols, col=ifelse(tttreat$logFC>0,"darkgreen","darkred"), las=3, cex.names = 1.2, cex.main = 1.2, axis.lty = 1, axes=T, ylab="logFC", ylim = c(-10,10))
abline(h=0)
dev.off()




