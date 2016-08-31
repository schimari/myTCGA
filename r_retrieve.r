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

# folder stuff
if(!dir.exists(download_folder)) { 
	dir.create(sprintf("%s", download_folder))
}
if(!dir.exists(sprintf("%s/%s", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s", download_folder, cdate))
}
if(!dir.exists(sprintf("%s/%s/data", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s/data", download_folder, cdate))
}

# rnaseqv2 data
if(!dir.exists(sprintf("%s/%s/data/rnaseq", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s/data/rnaseq", download_folder, cdate))
}

downloadTCGA(
	cancerTypes = ctype,
	dataSet = "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3",
	destDir = sprintf("%s/%s/data/rnaseq", download_folder, cdate),
	date = cdate
)

# clinical data
if(!dir.exists(sprintf("%s/%s/data/clinical", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s/data/clinical", download_folder, cdate))
}

downloadTCGA(
	cancerTypes = ctype,
	dataSet = "Merge_Clinical.Level_1",
	destDir = sprintf("%s/%s/data/clinical", download_folder, cdate),
	date = cdate
)

# oncotated mutation data
if(!dir.exists(sprintf("%s/%s/data/mutation", download_folder, cdate))) { 
	dir.create(sprintf("%s/%s/data/mutation", download_folder, cdate))
}

downloadTCGA(
	cancerTypes = ctype,
	dataSet = "Mutation_Packager_Oncotated_Calls.Level",
	destDir = sprintf("%s/%s/data/mutation", download_folder, cdate),
	date = cdate
)

# oncotated: 
#	tail -q -n +5 *.maf.txt >> ../oncotated.merged.maf.txt
#	copy and paste header line from z_oncotated.colnames.txt into file and save

# normal: merge maf files into one file with: cat *.maf.txt >> merged.maf.txt (e.g. from MINGW64 console)





	
