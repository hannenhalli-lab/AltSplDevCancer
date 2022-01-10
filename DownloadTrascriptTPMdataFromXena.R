args = commandArgs(trailingOnly=TRUE)
library(glue)
if (length(args) != 4) {
  stop("1. Provide TCGA CancerType ------ 2. Provide  TissueType ------ 3. logBase ------ 4. pseudoCount")
}


library(UCSCXenaTools)
library(dplyr)
library(TCGAbiolinks)
#library(preprocessCore)
#library(oncomix)
#source("/Users/singha30/DataAndAnalysis/barcodesFilterFunctionFromShixiangWang.R")



CancerType = args[1]
TissueType = args[2]
logBase = args[3]
pseudoCount = args[4]


antilogFunction <- function(x, b, p = pseudoCount) {
	alog <- NA
	tryCatch({
		alog <- (b ^ x) - p
	},error=function(x){})
	return (alog)
}																


CancerType <- cancerTypes[13]
outputFolder = CancerType
download <- "TPMvalues"
unit = "tpm"
if (nchar(CancerType) > 1) {
	if (file.exists(paste(CancerType, "/", sep = "/", collapse = "/"))) {
	    cat(paste(CancerType, "directory already exists...nothing created"))
	} else {
	    cat(paste(CancerType, "directory doesnt exists...creating new"))
	    dir.create(file.path(CancerType))
	}
	
	CancerType <- paste("TCGA-", CancerType, sep = "")
	query <- GDCquery(project = c(CancerType),
	               #data.category = "Sequencing Reads",
				   data.category = "Raw sequencing data", ####with legacy archive == T only
				   experimental.strategy = "RNA-Seq",
				   data.format = "FASTQ",  ####with legacy archive == T only
	               #sample.type = "Primary Tumor",
				   legacy = TRUE
				   )  
			   
	Samples <- getResults(query, cols = c("file_name","cases", "file_id", "sample_type", "is_ffpe"))
	Samples <- Samples[!Samples$is_ffpe, ]
	#FilteredSamples <- unique(Samples$cases)
	#########################clean the samples######################						
	#FilteredSamples <- tcga_replicateFilter(Samples$cases)
	#Samples <- Samples[Samples$cases %in% FilteredSamples, ]
	#Samples$cases <- gsub("\\w?-\\w+-\\w+-\\w+$", "", Samples$cases)
	
	SamplesToUse <- unique(gsub("\\w?-\\w+-\\w+-\\w+$", "", Samples$cases))
	
}		

subSets = 40
if (nchar(CancerType) > 1) {
	Selection = XenaData %>%
	  filter(XenaHostNames == "toilHub", 
	  		DataSubtype == "transcript expression RNAseq", 
			#XenaDatasets == "tcga_Kallisto_est_counts",
			XenaDatasets == glue("tcga_Kallisto_{unit}"),
			#XenaCohorts == "TARGET Pan-Cancer (PANCAN)"
			XenaCohorts == "TCGA Pan-Cancer (PANCAN)"
			#XenaDatasets == "tcga_rsem_isoform_tpm",
			)
		
	  Probes <- fetch_dataset_identifiers(Selection$XenaHosts, Selection$XenaDatasets)
	  
	  complete <- floor(length(Probes) / subSets)
	  
	  end = 0
	  cat("Fetching and writing file after taking antilog and subtracting pseudocount")	
	  for (i in 1:subSets) {
		  start = end + 1
		  end = end + complete + 1
		  FetchedExpressionValues <- fetch_dense_values(host = Selection$XenaHosts,
			                             dataset = Selection$XenaDatasets,
										 identifiers = Probes[start:end],
			  						  	 #samples = Samples$cases,
										 samples = SamplesToUse,
										 time_limit = 300,
			                             use_probeMap = F) #%>% data.frame()
	
		  sampleCount = ncol(FetchedExpressionValues)
		  naSums <- apply(FetchedExpressionValues,1,function(x) sum(is.na(x)))
		  indexes <- naSums == sampleCount
		  #if (indexes > 0) {
		  	FetchedExpressionValues <- FetchedExpressionValues[!indexes, ]
		  #}
		  if(i == 1){
			  write.table(FetchedExpressionValues, paste(outputFolder, "/Kallisto", download, "ForTranscriptExpressionIn", outputFolder, sep = ""), 
							sep = "\t", quote = F, row.names = T, col.names = T)
		  } else {
			  write.table(FetchedExpressionValues, paste(outputFolder, "/Kallisto", download, "ForTranscriptExpressionIn", outputFolder, sep = ""), 
							sep = "\t", quote = F, row.names = T, col.names = F, append = T)
		  }
		  message(paste("done for iteration no. ", i, sep = ""))
		  gc()
	  }
	  start = end + 1
	  end = length(Probes)
	  
	  if ((end - start) > 0 ) {
		  FetchedExpressionValues <- fetch_dense_values(host = Selection$XenaHosts,
						                             dataset = Selection$XenaDatasets,
													 identifiers = Probes[start:end],
						  						  	 samples = Samples$cases,
													 time_limit = 300,
						                             use_probeMap = F) #%>% collect()
	
		  
		  sampleCount = ncol(FetchedExpressionValues)
		  naSums <- apply(FetchedExpressionValues,1,function(x) sum(is.na(x)))
		  indexes <- naSums == sampleCount
		  #if (indexes > 0) {
		  	FetchedExpressionValues <- FetchedExpressionValues[!indexes, ]
			#}
		  write.table(FetchedExpressionValues, paste(outputFolder,"/Kallisto", download, "ForTranscriptExpressionIn", outputFolder, sep = ""), 
						sep = "\t", quote = F, row.names = T, col.names = F, append = T)
	  }	
}

if (nchar(TissueType) > 1) {
	
	if (file.exists(paste(TissueType, "/", sep = "/", collapse = "/"))) {
	    cat(paste(TissueType, "directory already exists...nothing created"))
	} else {
	    cat(paste(TissueType, "directory doesnt exists...creating new"))
	    dir.create(file.path(TissueType))
	}
	
	
	
	sampleAnnotations <- read.table("/Users/singha30/DataAndAnalysis/GtexSampleToTissueMap", sep = "\t", header = T)
	#sampleAnnotations <- read.table("GtexSampleToTissueMap", sep = "\t", header = T)
	samplesToConsider <- sampleAnnotations[sampleAnnotations$SMTSD == TissueType, 'SAMPID']
	
	Selection = XenaData %>% filter(XenaHostNames == "toilHub", DataSubtype == "transcript expression RNAseq", 
	  		XenaCohorts == "GTEX", XenaDatasets == glue("gtex_Kallisto_{unit}")
			)
	  
			
	  Probes <- fetch_dataset_identifiers(Selection$XenaHosts, Selection$XenaDatasets)
	  complete <- floor(length(Probes) / subSets)
	  
	  end = 0
	  cat("Fetching and writing file after taking antilog and subtracting pseudocount")	
	  for (i in 1:subSets) {
		  start = end + 1
		  end = end + complete + 1
		  FetchedExpressionValues <- fetch_dense_values(host = Selection$XenaHosts,
			                             dataset = Selection$XenaDatasets,
										 identifiers = Probes[start:end],
			  						  	 samples = as.character(samplesToConsider),
										 time_limit = 300,
			                             use_probeMap = F) #%>% data.frame()
	
		  
		  sampleCount = ncol(FetchedExpressionValues)
		  naSums <- apply(FetchedExpressionValues,1,function(x) sum(is.na(x)))
		  indexes <- naSums == sampleCount
		  #if (indexes > 0) {
		  	FetchedExpressionValues <- FetchedExpressionValues[!indexes, ]
		  #}
		  if(i == 1){
			   write.table(FetchedExpressionValues, paste(TissueType, "/Kallisto", download, "ForTranscriptExpressionIn", TissueType, sep = ""), 
							sep = "\t", quote = F, row.names = T, col.names = T)
		  } else {
			   write.table(FetchedExpressionValues, paste(TissueType, "/Kallisto", download, "ForTranscriptExpressionIn", TissueType, sep = ""), 
							sep = "\t", quote = F, row.names = T, col.names = F, append = T)
		  }
		  message(paste("done for iteration no. ", i, sep = ""))
		  gc()
	  }
	  start = end + 1
	  end = length(Probes)
	  
	  if ((end - start) > 0 ) {
		  FetchedExpressionValues <- fetch_dense_values(host = Selection$XenaHosts,
						                             dataset = Selection$XenaDatasets,
													 identifiers = Probes[start:end],
						  						  	 samples = Samples$cases,
													 time_limit = 300,
						                             use_probeMap = F) #%>% collect()
	
		  
		  sampleCount = ncol(FetchedExpressionValues)
		  naSums <- apply(FetchedExpressionValues,1,function(x) sum(is.na(x)))
		  indexes <- naSums == sampleCount
		  #if (indexes > 0) {
		  	FetchedExpressionValues <- FetchedExpressionValues[!indexes, ]
			#}
		  write.table(FetchedExpressionValues, paste(TissueType, "/Kallisto", download, "ForTranscriptExpressionIn", TissueType, sep = ""), 
						sep = "\t", quote = F, row.names = T, col.names = F, append = T)
	  }	
}
