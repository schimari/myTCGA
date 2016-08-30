# retrieve TCGA data

library("RTCGA")

# configuration
ctype <- "KIRC"
cdate <- tail(checkTCGA('Dates'), n=1)	# latest release
download_folder <- sprintf("%s", ctype)

# information
#infoTCGA()
#checkTCGA('Dates')
#checkTCGA('DataSets', cancerType = ctype)

# get data (rnaseqv2 and clinical)
if(!dir.exists(download_folder)) { 
	dir.create(sprintf("%s", download_folder))
}
if(!dir.exists(sprintf("%s/%s", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s", download_folder, cdate))
}
if(!dir.exists(sprintf("%s/%s/rnaseq", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s/rnaseq", download_folder, cdate))
}

downloadTCGA(
	cancerTypes = ctype,
	dataSet = "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3",
	destDir = sprintf("%s/%s/rnaseq", download_folder, cdate),
	date = cdate
)

if(!dir.exists(download_folder)) { 
	dir.create(sprintf("%s", download_folder))
}
if(!dir.exists(sprintf("%s/%s", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s", download_folder, cdate))
}
if(!dir.exists(sprintf("%s/%s/clinical", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s/clinical", download_folder, cdate))
}

downloadTCGA(
	cancerTypes = ctype,
	dataSet = "Merge_Clinical.Level_1",
	destDir = sprintf("%s/%s/clinical", download_folder, cdate),
	date = cdate
)

