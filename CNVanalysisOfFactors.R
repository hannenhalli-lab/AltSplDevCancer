library(TCGAbiolinks); library(glue); 
library(dplyr); library(data.table);
library(maftools); library(EnsDb.Hsapiens.v86)
library(GenomicFeatures); library(magrittr)
library(TCGAbiolinks); library(OneR)
EnsdB <- EnsDb.Hsapiens.v86
library(ggplot2); library(ggpubr)
#library(reshape2)


source('helperFunctions.R')
#cancerTypes = c("GBM", "LGG", "LIHC")

cnvResultList <- list()
Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),
GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))

preprocessCNVs <- function(CNVs, GenesDf) {
	CNVs <- CNVs %>% data.frame(check.names = F) %>% set_rownames(.[, 1]) %>% .[ ,-c(1:3)]
	rownames(CNVs) <- gsub("\\.\\d+$", "", rownames(CNVs))
	tumorBarcodes <- TCGAquery_SampleTypes(names(CNVs), typesample = c("TP"))
	CNVs <- CNVs[ ,colnames(CNVs) %in% tumorBarcodes]
	colnames(CNVs) <- gsub("\\-\\w+\\-\\w+\\-\\w+$", "", colnames(CNVs))
	
	geneNames <- GenesDf$GeneName[match(rownames(CNVs), GenesDf$GeneID)]
	cbind(geneNames, CNVs)
	#return(CNVs)
}


mapCNVs <- function(CNVs, sampleTypeDf, candidateFactor, returnType) {
	subData <- CNVs %>% .[CNVs$geneNames %in% candidateFactor, -1]
	names(subData) <- gsub("[A-Z]$", "", names(subData))
	gLevels <- sampleTypeDf$sampleType %>% unique %>% as.character()
	gainList <- list()
	for (index in 1:length(gLevels)) {
		gLevel <- gLevels[index]
		samples <- sampleTypeDf %>% subset(., sampleType == gLevel, select = sample) %>% unlist() %>% as.character()
		if (returnType == "gaincount") {
			gain <- subData %>% .[ ,names(.) %in% samples] %>% apply(., 2, function(x) sum(x == 1))
			gainList[[index]] <- gain
		}
		if (returnType == "netgain") {
			netGain <- subData %>% .[ ,names(.) %in% samples] %>% apply(., 2, function(x) sum(x))
			gainList[[index]] <- netGain
		}
	}
	gainList <- unlist(gainList)
	sampleTypeDf <- merge(sampleTypeDf, gainList, by.x = "sample", by.y = 0)
	if (returnType == "gaincount") {
		names(sampleTypeDf)[names(sampleTypeDf) == "y"] <- "gaincount"
	} else {
		names(sampleTypeDf)[names(sampleTypeDf) == "y"] <- "netgain"
	}
	return(sampleTypeDf)
}


getEvents <- function (AssociationClusters, ExonType, Association) {
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	ExonType <- paste(ExonType, Association, sep = ".")
	events <- AssociationClusters[AssociationClusters$ExonType == ExonType, 'Exon'] %>% as.character()
	return(events)
}



factorWiseAnalysis <- function(inputFile, candidateFactors, cancerGeneTPM, cnvType, binFactor, posPSI) {
	outputList <- list()
	for (factorIndex in 1:length(candidateFactors)) {
		candidateFactor <- candidateFactors[factorIndex]
		medianExpression <- cancerGeneTPM %>% .[rownames(.) %in% candidateFactor, ] %>% unlist()#%>% apply(., 2, median, na.rm = T, na.action = na.pass)
		nbins <- 3; binNames <- paste("Level", 1:nbins, sep = "")
		#expressionBins <- bin(medianExpression, nbins, labels = binNames, method = "length")
		
		if (binFactor == "CNV") {
			subData <- CNVs %>% .[CNVs$geneNames %in% candidateFactor, -1] %>% unlist()
			names(subData) <- gsub("[A-Z]$", "", names(subData))
			#sampleBins = c(-1, 0, 1)
			#sampleBins = sapply(1:3, function(index) {subData == sampleBins[index]})
			sampleTypeDf <- data.frame(sample = names(subData), sampleType = subData)
			sampleTypeDf <- merge(sampleTypeDf, posPSI, by.x = "sample", by.y = 0)
			sampleTypeDf <- merge(sampleTypeDf, medianExpression, by.x = "sample", by.y = 0)
			names(sampleTypeDf)[names(sampleTypeDf) == "y.x"] <- "inclusionLevel"
			names(sampleTypeDf)[names(sampleTypeDf) == "y.y"] <- "expression"
			
			mappedCNVs <- sampleTypeDf %>% group_by(sampleType) %>% summarise_at(vars("inclusionLevel", "expression"), list(mean)) %>% data.frame()
			outputList[[candidateFactor]] <- mappedCNVs
		}else {
			if (binFactor == "splicing") {
				#nbins <- 2; binNames <- paste("Level", 1:nbins, sep = "")
				#sampleBins <- bin(posPSI, nbins, labels = binNames, method = "content")
				#sampleBins <- bin(posPSI, nbins, labels = binNames, method = "cluster")
				
				level1 = summary(posPSI)['1st Qu.']
				level2 = summary(posPSI)['3rd Qu.']
				level3 = summary(posPSI)['Max.']
				     
				sampleBins <- cut(posPSI, c(0, level1, level2, level3), labels = c("Level1", "Level2", "level3"))
				#sampleBins <- cut(posPSI, c(0,0.30,1), labels = c("Level1", "Level2"))
			
				sampleTypeDf <- data.frame(sample = names(posPSI), inclusionLevel = posPSI, sampleType = sampleBins)
				sampleTypeDf <- merge(sampleTypeDf, medianExpression, by.x = "sample", by.y = 0)
				names(sampleTypeDf)[names(sampleTypeDf) == "y"] <- "expression"
			}
			if (binFactor == "expression") {
				sampleBins <- NA; mappedCNVs <- NA
				tryCatch({
					sampleBins <- bin(medianExpression, nbins, labels = binNames, method = "content")
			    },error=function(x){})
				
				if (length(sampleBins) > 1)	{
					sampleTypeDf <- data.frame(sample = names(medianExpression), expression = medianExpression, sampleType = sampleBins)
					sampleTypeDf <- merge(sampleTypeDf, posPSI, by.x = "sample", by.y = 0)
					names(sampleTypeDf)[names(sampleTypeDf) == "y"] <- "inclusionLevel"
				}else {
					
				}
			}
			mappedCNVs <- mapCNVs(CNVs = inputFile, sampleTypeDf = sampleTypeDf, candidateFactor = candidateFactor, returnType = cnvType) %>% 
							group_by(sampleType) %>% summarise_at(vars(glue(cnvType), "inclusionLevel", "expression"), list(mean)) %>% data.frame()
					
			outputList[[candidateFactor]] <- mappedCNVs
		}
		
	}
	return(outputList)
}


cancerType <- "LIHC"; cancerType <- "LIHC"; devTissue <- "Liver"

load(glue('plsrModelForMedianPosEmbSplicingIn{devTissue}.rda')) ####model for splicing factors
load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda")) #### data for splicing events
embPos <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Positive")



pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/Cancer/Kallisto/"
cancerExpression <- glue("{pathToExpressionFile}TCGA/Expression/KallistoTPMvaluesForTranscriptExpressionIn{cancerType}")
cancerPSIfile <- glue("{pathToExpressionFile}TCGA/PSIvalues/SuppaPSIvaluesUsingKallistoIn{cancerType}.psi")

cancerExpression <- read.table(cancerExpression, sep = "\t", header = T, check.names = F)
names(cancerExpression) <- gsub("\\.", "-", names(cancerExpression))
cancerSamples <- TCGAquery_SampleTypes(names(cancerExpression), typesample = c("TP"))
cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples]
cancerGeneTPM <- processEpxressionFile(dataFrame = cancerExpression, Transcripts = Transcripts, replicates = F)
rownames(cancerGeneTPM) <- gsub("-", "\\.", rownames(cancerGeneTPM))

cancerPSI <- read.table(cancerPSIfile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)
names(cancerPSI) <- gsub("\\.", "-", names(cancerPSI))
cancerPSI <- naFilter(dataFrame = cancerPSI, cutoff = 0.5)
cancerPSI <- cancerPSI[, names(cancerPSI) %in% cancerSamples]

#embNeg <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Negative", 'Exon'] %>% as.character()
posPSI <- cancerPSI[rownames(cancerPSI) %in% embPos, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)

	
	
##########CNV analysis##########
load(glue('~/DataAndAnalysis/TissueRelevance/CancerWiseCNVs/{cancerType}gisticCNVscores.rda'))
CNVs <- preprocessCNVs(CNVs = data, GenesDf = GenesDf)



#candidateFactors <- factorsDf %>% rownames(.) %>% .[1:100]
candidateFactors <- combinedLoadings %>% rownames(.)
criticalFactors <- combinedLoadings %>% .[.$FDR <= 0.05 & .$coefficient > 0, ] %>% .[order(.$coefficients, decreasing = T), ] %>% rownames(.) %>% .[1:100]
nonCriticalFactors <- combinedLoadings %>% .[.$coefficients <= 0.0, ] %>% rownames()



#cnvType = "gaincount"
cnvType = "netgain"
binFactor = "expression"
#binFactor = "splicing"


copyNumberAndInclusion <- factorWiseAnalysis(inputFile = CNVs, candidateFactors = candidateFactors, cancerGeneTPM = cancerGeneTPM, cnvType = cnvType, binFactor = binFactor, posPSI = posPSI)
copyNumberAndInclusion <- do.call("rbind", copyNumberAndInclusion)
copyNumberAndInclusion$GeneName <- gsub("\\.\\d+$", "", rownames(copyNumberAndInclusion))
copyNumberAndInclusion$GeneName <- factor(copyNumberAndInclusion$GeneName, levels = unique(copyNumberAndInclusion$GeneName))
#ggplot(data = copyNumberAndInclusion, aes(x = sampleType, y = netgain)) + geom_bar(stat = "identity") + facet_wrap(~GeneName, scales = "free_y")

copyNumberAndInclusion$factorType <- "Other"
copyNumberAndInclusion$factorType[copyNumberAndInclusion$GeneName %in% criticalFactors] = "Critical"
copyNumberAndInclusion$factorType[copyNumberAndInclusion$GeneName %in% nonCriticalFactors] = "Non critical"
copyNumberAndInclusion <- copyNumberAndInclusion[copyNumberAndInclusion$factorType!= "Other", ]


cnvResultList[[devTissue]] <- copyNumberAndInclusion

save(cnvResultList, file = "cnvResultsForCriticalFactorsInCancerTypes.Rda")


