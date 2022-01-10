library(dplyr); library(magrittr); library(ensembldb); 
library(EnsDb.Hsapiens.v86); library(clusterProfiler)
library(gtools); library(reshape2); library(ggplot2)
library(clusterProfiler); library(gplots); library(glue)
orgDb <- "org.Hs.eg.db"; EnsdB <- EnsDb.Hsapiens.v86
library('rstatix'); library(ggpubr)
library(data.table)
source('../helperFunctions.R')
library(pls) ### for using predict function

Transcripts <- transcripts(EnsdB,
			          columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "tx_seq_start", "tx_seq_end"), #listColumns(EnsdB , "tx")),
			          #filter = TxBiotypeFilter("nonsense_mediated_decay"),
			          return.type = "DataFrame")
GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name))


getEvents <- function (AssociationClusters, ExonType, Association) {
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	ExonType <- paste(ExonType, Association, sep = ".")
	events <- AssociationClusters[AssociationClusters$ExonType == ExonType, 'Exon'] %>% as.character()
	return(events)
}


sampleInfo <- read.csv("sample_info.csv", header = T)
annotations <- fread("Cell_lines_annotations_20181226.txt", sep = "\t", header = T) %>% data.frame() %>% 
	mutate(depMapID = gsub("-", "\\.", depMapID)) %>% subset(. , select = c(depMapID, Doubling.Time.from.Vendor, Doubling.Time.Calculated.hrs))





devTissues <- c("Hindbrain" = "Brain Cancer", Kidney = "Kidney Cancer", Liver = "Liver Cancer")
desiredTypes <- list(Hindbrain = "all", Kidney = "all", Liver = "all")

dataForplotList <- list()
for (devTissue in names(devTissues)) {
	
	cancer <- devTissues[devTissue]
	cancerType <- gsub(" ", "", cancer)
	desiredType <- desiredTypes[[devTissue]]
	
	load(glue("../../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))
	positiveExons <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Positive") %>% gsub(".\\d+;SE", ";SE", .)
	negativeExons <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Negative") %>% gsub(".\\d+;SE", ";SE", .)
	
	expressionFile <- glue("ccleRSEMtranscriptTPMfor{cancerType}")
	psiFile <- glue("SuppaPSIvaluesIn{cancerType}Lines.psi")


	expressionFile <- read.table(expressionFile, sep = "\t", header = T)
	psiFile <- read.table(psiFile, sep = "\t", header = T, na.strings = c('nan', 'NA'))

	geneLevelTPM <- processEpxressionFile(dataFrame = expressionFile, Transcripts = Transcripts, replicates = F)
	psiFile <- naFilter(psiFile, cutoff = 0.5)


	medianPosSplicing <- psiFile[rownames(psiFile) %in% positiveExons, ] %>% apply(., 2, median, na.rm = T, na.action = na.pass)
	medianNegSplicing <- psiFile[rownames(psiFile) %in% negativeExons, ] %>% apply(., 2, median, na.rm = T, na.action = na.pass)

	
	samplesToUse <- names(geneLevelTPM) 

	subSamples <- sampleInfo[, c('DepMap_ID', 'Subtype', 'primary_or_metastasis', 'primary_disease')] %>% 
			mutate(DepMap_ID = gsub("-", "\\.", DepMap_ID)) %>% subset(., DepMap_ID %in% samplesToUse)
	
	subSamples <- merge(subSamples, annotations, by.x = "DepMap_ID", by.y = "depMapID")
	
		
	if (desiredType == "all") {
		subSamples <- subSamples[subSamples$primary_disease == cancer, ]
	} else {
		subSamples <- subSamples[subSamples$primary_disease == cancer, ] %>% .[.$Subtype %in% desiredType, ]
	}

	load(glue("../plsrModelForMedianPosEmbSplicingIn{devTissue}.rda"))

	rownames(geneLevelTPM) <- gsub("-", "\\.", rownames(geneLevelTPM))
	testMatrix <- geneLevelTPM[rownames(geneLevelTPM) %in% rownames(combinedLoadings), ] %>% t() %>% data.frame()


	testMatrix <- log2(testMatrix + 1) %>% as.matrix()
	testVariables <- data.frame(posPSI = medianPosSplicing) %>% as.matrix()
	#testVariables <- scale(testVariables, center = T, scale = F)
	predictedValues <- predict(model, ncomp = 2, newdata = testMatrix)
	
	predictedValuesDf <- data.frame(cancer = cancer, actualValues = testVariables, predictedValues = predictedValues)
	dataForplotList[[cancer]] <- predictedValuesDf
}
save(dataForplotList, file = "actualAndpredictedEmbryonicPositiveSplicingInCCLE.Rda")



