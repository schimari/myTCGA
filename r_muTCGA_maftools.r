# maftools: https://github.com/PoisonAlien/maftools
# data from: http://firebrowse.org/?cohort=KIRC&download_dialog=true
# or from: https://gdc-docs.nci.nih.gov/Data/Release_Notes/Data_Release_Notes/	<- mutec somatic
# or Oncotated Calls with protein change: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRC/20160128/gdac.broadinstitute.org_KIRC.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0.tar.gz

library("maftools")
library("ggrepel")

# configure
cohort <- "KIRC"
fdate <- "2016-01-28"
bfolder <- sprintf("%s/%s", cohort, fdate)

# oncotated data
data.maf <- read.maf(maf = sprintf("%s/data/mutation/oncotated.merged.maf.txt", bfolder), removeSilent = T, useAll = F, isTCGA = T)

# output and plot folder
if(!dir.exists(sprintf("%s/maftools_output", bfolder))) { 
	dir.create(sprintf("%s/maftools_output", bfolder))
}
if(!dir.exists(sprintf("%s/maftools_plots", bfolder))) { 
	dir.create(sprintf("%s/maftools_plots", bfolder))
}

# MAF file summarized by Genes 
getGeneSummary(data.maf)

# MFA file summarized by samples (Tumor Sample Barcode)
getSampleSummary(data.maf)

# Writes maf summary to an output file with basename data.maf.
write.mafSummary(maf = data.maf, basename = sprintf("%s/maftools_output/%s_maf", bfolder, cohort))

# Quicky plot maf stats
plotmafSummary(maf = data.maf, rmOutlier = T, addStat = 'median', dashboard = TRUE, 
	file=sprintf("%s/maftools_plots/%s_mafSummary.pdf", bfolder, cohort), width = 11.69, height = 8.27)

# oncoplot to summarize maf file
cairo_pdf(sprintf("%s/maftools_plots/%s_oncoplot.pdf", bfolder, cohort), width=11.69, height=8.27)
	oncoplot(maf = data.maf, top = 10, removeNonMutated = TRUE)
dev.off()

# Classify SNVs into Trasitions and Transversions
cairo_pdf(sprintf("%s/maftools_plots/%s_titv.pdf", bfolder, cohort), width=11.69, height=8.27)
	data.titv = titv(maf = data.maf, useSyn = T)
dev.off()

# Lollipop plots for amino acid changes	-> only works with oncotated calls
int_gene <- "TTN"
cairo_pdf(sprintf("%s/maftools_plots/%s_lollipop_%s.pdf", bfolder, cohort, int_gene), width=11.69, height=8.27)
	lollipopPlot(maf = data.maf, gene = int_gene, AACol = "Protein_Change", labelPos = "all")
dev.off()

# Genecloud
cairo_pdf(sprintf("%s/maftools_plots/%s_genecloud.pdf", bfolder, cohort), width=11.69, height=8.27)
	geneCloud(input = data.maf, minMut = 10, top = 20)
dev.off()
	
# Detcting cancer causing genes	-> only works with oncotated calls
data.sig = oncodrive(maf = data.maf, AACol = "Protein_Change", minMut = 5, pvalMethod = "zscore")
cairo_pdf(sprintf("%s/maftools_plots/%s_oncodrive.pdf", bfolder, cohort))
	plotOncodrive(res = data.sig, fdrCutOff = 0.1, useFraction = T)
dev.off()




