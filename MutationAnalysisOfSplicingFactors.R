library(TCGAbiolinks); library(glue); 
library(dplyr); library(data.table);
library(maftools); library(EnsDb.Hsapiens.v86)
library(GenomicFeatures); library(magrittr)
library(TCGAbiolinks); library(OneR)
EnsdB <- EnsDb.Hsapiens.v86
library(ggplot2); library(ggpubr)
library(scales)
#library(reshape2)


source('helperFunctions.R')
#CancerTypes = c("GBM", "LGG", "LIHC")


Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),
GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))


#transcriptLengthVec <- transcriptLengths(EnsdB)
#transcriptLengthVec <- transcriptLengthVec %>% group_by(gene_id) %>% summarise(., length = max(tx_len))
#GenesDf <- merge(GenesDf, transcriptLengthVec, by.x = "GeneID", by.y = "gene_id")


#GenesDf <- genes(EnsdB, GeneBiotypeFilter("protein_coding"), columns = c("gene_id", "seq_name", "gene_name", "seq_strand"), return.type = "data.frame")

getMutationData <- function(CancerType, pipelines, download, filePath) {
	if (download) {
		query <- 1; ntries = 0
		maf <- NA
		while(length(maf) == 1 & ntries < 10) {
			tryCatch({
				#query <- GDCquery(project = project, data.category = dataCategory, data.type = dataType, workflow.type = workFlow)
				maf <- GDCquery_Maf(CancerType, pipelines = "mutect2")
			    },error=function(x){})
			ntries <- ntries + 1
			if (ntries == 10) {
				print(glue("Could not fectch data even after 10 attempts for {CancerType}"))
			}
		}
		fwrite(maf, filePath, quote = "auto", sep = "\t")
	}else {
		maf <- read.csv(filePath, sep = "\t", header = T)
	}
	return(maf)
}



getEvents <- function (AssociationClusters, ExonType, Association) {
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	ExonType <- paste(ExonType, Association, sep = ".")
	events <- AssociationClusters[AssociationClusters$ExonType == ExonType, 'Exon'] %>% as.character()
	return(events)
}




#cancerType <- c("LGG", "GBM"); devTissue <- "Hindbrain"
#cancerType <- c("LGG"); devTissue <- "Hindbrain"
#cancerType <- c("LIHC"); devTissue <- "Liver"

devTissue <- c("Hindbrain", "Liver", "Kidney")
cancerType <- c("Hindbrain" = "LGG", "Liver" = "LIHC", "Kidney" = "KIRP")

download = F

mafFiles <- list()
if (download) {
	for (cancer in cancerType) {
		filePath = glue("MAFs/{cancer}mutationsFromTCGA.maf")
		mafFiles[[cancer]] <- getMutationData(CancerType = cancer, pipelines = "mutect2", download = T, filePath = filePath)
	}
	
} else {
	for (cancer in cancerType) {
		filePath = glue("MAFs/{cancer}mutationsFromTCGA.maf")
		mafFiles[[cancer]] <- read.maf(filePath)
	}
}
mutationDataList <- list(); maxLoad <- 500

for (cancer in cancerType) {
	mutationSubData <- mafFiles[[cancer]]@data %>% subset(., Consequence != "synonymous_variant" & Mutation_Status == "Somatic", select = c(Hugo_Symbol, Mutation_Status, Variant_Classification, Consequence, Tumor_Sample_Barcode, SIFT, PolyPhen))
	mutationSubData$SIFT <- gsub("\\(.+", "", mutationSubData$SIFT)
	mutationSubData$PolyPhen <- gsub("\\(.+", "", mutationSubData$PolyPhen)
	#mutationDf <- mutationSubData[grep("deleterious", mutationSubData$SIFT), ]
	hyperMutator <- mutationSubData %>% group_by(., Tumor_Sample_Barcode) %>% summarize(mutationLoad = length(Tumor_Sample_Barcode))
	desiredSamples <- hyperMutator %>% subset(., mutationLoad < maxLoad, select = Tumor_Sample_Barcode) %>% unlist
	mutationSubData <- mutationSubData[mutationSubData$Tumor_Sample_Barcode %in% desiredSamples, ]
	
	mutationDataList[[cancer]] <- mutationSubData
}


###description of cnv data####
##https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
cnvDataList <- list()
for (cancer in cancerType) {
	load(glue('~/DataAndAnalysis/TissueRelevance/CancerWiseCNVs/{cancer}gisticCNVscores.rda'))
	cnvDataList[[cancer]] <- preprocessCNVs(data, GenesDf)
}



#mutationSubData <- do.call("rbind", mutationDataList)
#rm(mutationDataList)


readoutList <- list()
regressionCoefficientList <- list()
tfList <- list()
for (index in 1:length(cancerType)) {
	cancer <- cancerType[index]
	tissue <- names(cancerType)[index]
	load(glue('plsrModelForMedianPosEmbSplicingIn{tissue}.rda')) ####model for splicing factors
	
	load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{tissue}.Rda")) #### data for splicing events
	embPos <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Positive")
	#embNeg <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Negative")
	
	tempDf <- combinedLoadings %>% .[sort(rownames(.)), ] %>% subset(., select = c(coefficients, FDR)) %>% mutate(factorName = rownames(.), tissue = tissue, cancer = cancer, .before = 1)
	regressionCoefficientList[[cancer]] <- tempDf
	readoutList[[cancer]] <- embPos
	
	if (file.exists(glue("TFbindingEnrichmentIn1kbPromotersUsingTFEA.ChIPremapDataIn{tissue}SplicingFactors.Rda"))) {
		load(glue("TFbindingEnrichmentIn1kbPromotersUsingTFEA.ChIPremapDataIn{tissue}SplicingFactors.Rda"))
		tfList[[cancer]] <- rankedTFs[1:20, 'TF']
	}
}


coeffcientsDf <- lapply(regressionCoefficientList, function(x) x$coefficients %>% set_names(x$factorName)) %>% do.call("cbind", .) 
#coeffcientsDf <- coeffcientsDf %>% normalize.quantiles() %>% data.frame() %>% set_names(cancerType) %>% set_rownames(rownames(coeffcientsDf)) 
fdrDf <- lapply(regressionCoefficientList, function(x) x$FDR %>% set_names(x$factorName)) %>% do.call("cbind", .) %>% data.frame() %>% set_names(cancerType)

cancerExpressionList <- list()
pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/Cancer/Kallisto/"
for (cancer in cancerType) {
	cancerExpression <- glue("{pathToExpressionFile}TCGA/Expression/KallistoTPMvaluesForTranscriptExpressionIn{cancer}")
	cancerExpression <- read.table(cancerExpression, sep = "\t", header = T, check.names = F)
	names(cancerExpression) <- gsub("\\.", "-", names(cancerExpression))

	cancerSamples <- TCGAquery_SampleTypes(names(cancerExpression), typesample = c("TP"))
	cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples]
	cancerExpression <- cancerExpression[order(rownames(cancerExpression)), ]
	#cancerExpressionList[[cancer]] <- cancerExpression
	
	cancerGeneTPM <- processEpxressionFile(dataFrame = cancerExpression, Transcripts = Transcripts, replicates = F)
	colnames(cancerGeneTPM) <- gsub("\\w+\\.", "", names(cancerGeneTPM))
	cancerExpressionList[[cancer]] <- cancerGeneTPM
	
}

#cancerExpression <- do.call("cbind", cancerExpressionList)
#cancerGeneTPM <- processEpxressionFile(dataFrame = cancerExpression, Transcripts = Transcripts, replicates = F)
#colnames(cancerGeneTPM) <- gsub("\\w+\\.", "", names(cancerGeneTPM))
#rm(cancerExpressionList)



cancerPSIlist <- list()
for (cancer in cancerType) {
	cancerPSIfile <- glue("{pathToExpressionFile}TCGA/PSIvalues/SuppaPSIvaluesUsingKallistoIn{cancer}.psi")
	cancerPSI <- read.table(cancerPSIfile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)
	names(cancerPSI) <- gsub("\\.", "-", names(cancerPSI))

	cancerPSI <- naFilter(dataFrame = cancerPSI, cutoff = 0.5)
	
	cancerSamples <- TCGAquery_SampleTypes(names(cancerPSI), typesample = c("TP"))
	cancerPSI <- cancerPSI[, names(cancerPSI) %in% cancerSamples]
	
	cancerPSI <- cancerPSI[order(rownames(cancerPSI)), ]
	cancerPSIlist[[cancer]] <- cancerPSI
}

#cancerPSI <- do.call("cbind", cancerPSIlist)
#rm(cancerPSIlist)

devTissues <- c("Hindbrain", "Liver", "Kidney")
cancerTypes <- c("Hindbrain" = "LGG", "Liver" = "LIHC", "Kidney" = "KIRP")


#intactExpression <- criticalExpression[names(criticalExpression) %in% intactSamples]


weightedDistance <- function(A, B, W) {
	coefSums <- sqrt(sum(W * (A - B)^2))
	return(coefSums)
}

#W = combinedLoadings$coefficients
#W = combinedLoadings$coefficients * ifelse(combinedLoadings$FDR < 0.05, 1, 0)


euclidDistanceFunction <- function(geneExpression, geneList, intactSamples, weightVector) {
	predictors <- geneExpression %>% .[rownames(.) %in% geneList, ] %>% .[order(rownames(.)), ]
	weightVector <- weightVector[order(weightVector)]
	
	if (length(weightVector) == length(geneList)) {
		vectorNames <- combn(colnames(predictors), 2)
		euclidDistances <- apply(vectorNames, 2, function(x) {A = x[1]; B = x[2]; A = predictors[ ,A]; B = predictors[ ,B]; W = weightVector; weightedDistance(A, B, W)})
		weightedDistances1 <- data.frame(Var1 = vectorNames[1, ], Var2 = vectorNames[2, ], value = euclidDistances)
		weightedDistances2 <- data.frame(Var1 = vectorNames[2, ], Var2 = vectorNames[1, ], value = euclidDistances)
		euclidDistances <- rbind(weightedDistances1, weightedDistances2)
	}
	
	if (weightVector == 1) {
		euclidDistances <- dist(t(predictors), method = "euclidean", diag = T) %>% as.matrix()
		euclidDistances[lower.tri(euclidDistances)] <- t(euclidDistances)[lower.tri(euclidDistances)]

		euclidDistances <- reshape2::melt(euclidDistances)
		euclidDistances <- euclidDistances[euclidDistances$Var1 != euclidDistances$Var2, ]
	}
	
	if (intactSamples == "All") {
		return(euclidDistances)
	} else {
		euclidDistances <- euclidDistances[euclidDistances$Var2 %in% intactSamples, ]
		return(euclidDistances)
	}
}


load("commonEventsAndFactorsInLiverBrainKidney.Rda")

#mutationClasses <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins")
mutationClasses <- c("Nonsense_Mutation")
mutationConsequence <- "NMD_|stop"
mutatedPSIdf <- NULL
mutationType <- "SNP"
useSpecific <- T; weighted = F
exonWise <- F
for (cancer in cancerType) {
	
	cancerGeneTPM <- cancerExpressionList[[cancer]]
	rownames(cancerGeneTPM) <- gsub("-", "\\.", rownames(cancerGeneTPM))
	
	
	if (mutationType == "SNP") {
		#mutationDf <- mutationDataList[[cancer]] %>% .[.$Variant_Classification == "Nonsense_Mutation" | .$SIFT == "deleterious", ]
		#mutationDf <- mutationDataList[[cancer]] %>% .[.$Variant_Classification == "Nonsense_Mutation", ]
		
		mutationDf1 <- mutationDataList[[cancer]] %>% .[.$Variant_Classification %in% mutationClasses, ]
		mutationDf2 <- mutationDataList[[cancer]] %>% .[grep(mutationConsequence, .$Consequence), ]
	
		mutationDf <- rbind(mutationDf1, mutationDf2) %>% unique()
	
		mutationDf$Tumor_Sample_Barcode <- gsub("\\-\\w+\\-\\w+\\-\\w+$", "", mutationDf$Tumor_Sample_Barcode)
		mutationDf$Tumor_Sample_Barcode <- gsub("[A-Z]$", "", mutationDf$Tumor_Sample_Barcode)
		mutationDf <- reshape2::dcast(mutationDf, Hugo_Symbol ~ Tumor_Sample_Barcode, fun.aggregate = length)
	
		mutationDf <- mutationDf %>% subset(., Hugo_Symbol %in% allFactors)
	
		mutatedSamples <- mutationDf[ ,-1] %>% colSums %>% .[. > 0] %>% names()
		intactSamples <- colnames(cancerGeneTPM)[!(colnames(cancerGeneTPM) %in% mutatedSamples)]
	} 
	if (mutationType == "CNV") {
		mutationDf <- cnvDataList[[cancer]] %>% subset(., geneNames %in% allFactors) %>% set_rownames(.$geneNames) %>% set_colnames(gsub("\\w$", "", colnames(.))) %>% na.omit()
		mutatedSamples <- mutationDf[ ,-1] %>% colSums %>% .[. != 0] %>% names()
		#intactSamples <- colnames(cancerGeneTPM)[!(colnames(cancerGeneTPM) %in% mutatedSamples)]
		intactSamples <- colnames(cancerGeneTPM) #### use all samples as intact samples for cnvs
	}
	
	########CHECK THIS for CNV #####
	#intactSamples <- intactSamples[!(intactSamples %in% cnvSamples)]
	
	if (weighted) {
		W = -coeffcientsDf %>% .[ , cancer] #### minus sign is used so that highly predictive factors may result in lower euclidian distance
		W <- rescale(W, to = c(0, 1)) %>% set_names(rownames(coeffcientsDf))
	} else {
		W <- 1
	}
	
	vectorDistDf <- euclidDistanceFunction(geneExpression = cancerGeneTPM, geneList = allFactors, intactSamples = intactSamples, weightVector = W)
	
	if (exonWise) {
		if (useSpecific) {
			events <- readoutList[[cancer]]
			readout <- cancerPSIlist[[cancer]] %>% .[rownames(.) %in% events, ] 
		}else {
			readout <- cancerPSIlist[[cancer]] %>% .[rownames(.) %in% commonEvents, ]
		}
	}else {
		if (useSpecific) {
			events <- readoutList[[cancer]]
			readout <- cancerPSIlist[[cancer]] %>% .[rownames(.) %in% events, ] %>% apply(., 2, median, na.rm = T, na.action = na.pass)
		}else {
			readout <- cancerPSIlist[[cancer]] %>% .[rownames(.) %in% commonEvents, ] %>% apply(., 2, median, na.rm = T, na.action = na.pass)
		}
	}
	
	
	
	for (index in 1:nrow(mutationDf)) {
		
		if (mutationType == "SNP") {
			factor = mutationDf[index, 'Hugo_Symbol']
			#factorType = mutationDf[index, 'factorType']
			tempDf <- mutationDf[ ,-(1:2)] %>% data.matrix
			mutations <- tempDf[index, ] %>% unlist
	
			#intactSamples <- mutations[mutations == 0] %>% names
			mutatedSamples <- mutations[mutations > 0] %>% names
			proceed <- T
		}
		if (mutationType == "CNV") {
			factor = mutationDf[index, 'geneName'] %>% as.character()
			
			#factorType = mutationDf[index, 'factorType']
			tempDf <- mutationDf[ ,-(1)] %>% data.matrix
			mutations <- tempDf[index, ] %>% unlist
	
			
			mutatedSamples <- mutations[mutations < 0] %>% names
			intactSamples <- mutations[mutations == 0] %>% names
			
			mutateExpression <- cancerGeneTPM %>% .[rownames(.) %in% factor , colnames(.) %in% mutatedSamples] %>% unlist() %>% median(., na.rm = T)
			intactExpression <- cancerGeneTPM %>% .[rownames(.) %in% factor , colnames(.) %in% intactSamples] %>% unlist() %>% median(., na.rm = T)
			
			if (intactExpression > 0) {
				foldChange <- mutateExpression / intactExpression
				print(foldChange)
			}else {
				foldChange <- 10
			}
			
			if (sum(foldChange < 0.5 & length(mutateExpression) > 0) == 1) {
				proceed <- T
			} else {
				proceed <- F
			}
		}		
		
		if (proceed) {
			regressionCoefficient <- coeffcientsDf %>% .[rownames(.) %in%  factor, cancer]
			fdr <- fdrDf %>% .[rownames(.) %in%  factor, cancer]
		
	
			for (mutatedSample in mutatedSamples) {
				subData <- vectorDistDf[vectorDistDf$Var1 %in% mutatedSample, ] %>% .[order(.$value, decreasing = F), ]
	
				similarSamples <- subData$Var2[1:10] %>% as.character()
				maxDistvalue <- max(subData$value[1:10])
	
				#similarExpression <- criticalExpression[names(criticalExpression) %in% similarSamples]
				#readoutSd <- sd(readout)##
				
				if (exonWise) {
					mutatedPSI <- NA; deltaPSImedians <- NA; medianMutatedPSI <- NA; 
					increaseUniverse <- NA; decreaseUniverse <- NA
					
					effectSize = 0.1
					
					expectedPSI <- readout[ ,names(readout) %in% similarSamples] %>% apply(., 1, median, na.rm = T, na.action = na.pass)
					mutatedPSI <- readout[ ,names(readout) == mutatedSample]
					deltaPSI <- mutatedPSI - expectedPSI
					
					medianExpectedPSI <- median(expectedPSI, na.rm = T)
					if (length(mutatedPSI) != 0) {
						medianMutatedPSI <- median(mutatedPSI, na.rm = T)
						deltaPSImedians <- medianMutatedPSI - medianExpectedPSI
						increaseUniverse <- expectedPSI %>% na.omit() %>% .[. < 0.9] %>% length
						decreaseUniverse <- expectedPSI %>% na.omit() %>% .[. > 0.1] %>% length
					}
					
					
					fractionIncreased <- sum(deltaPSI > effectSize, na.rm = T) / increaseUniverse
					fractionDecreased <- sum(deltaPSI < -effectSize, na.rm = T) / decreaseUniverse
					
					foldChange <- log2(fractionIncreased / fractionDecreased)
					
					#effectClass <- ifelse(foldChange)
					
					tempDf <- c("factor" = factor, "cancer" = cancer, "regressionCoefficient" = regressionCoefficient, "fdr" = fdr, "maxDistValue" = maxDistvalue, 
								"medianExpectedPSI" = medianExpectedPSI, medianMutatedPSI = medianMutatedPSI, "deltaOfMedians" = deltaPSImedians, 
								"fractionIncreased" = fractionIncreased, "fractionDecreased" = fractionDecreased, 
								"increaseUniverse" = increaseUniverse, "decreaseUniverse" = decreaseUniverse , "foldChange" = foldChange)
								
					mutatedPSIdf <- rbind(mutatedPSIdf, tempDf)
					
				} else {
					expectedPSI <- readout[names(readout) %in% similarSamples] %>% median(., na.rm = T, na.action = na.pass)
					minExp <- readout[names(readout) %in% similarSamples] %>% min(., na.rm = T)
					maxExp <- readout[names(readout) %in% similarSamples] %>% max(., na.rm = T)
				
				
					deviationPSI <- readout[names(readout) %in% similarSamples] %>% sd(., na.rm = T)
					mutatedPSI <- readout[names(readout) == mutatedSample]%>% median(., na.rm = T, na.action = na.pass)
	
					if (is.na(mutatedPSI)) {
						#tempDf <- data.frame("factor" = factor, "factorType" = factorType, "expectedPSI" = expectedPSI, "mutatedPSI" = NA, row.names = NULL)
					}else {
						#tempDf <- data.frame("factor" = factor, "factorType" = factorType, "expectedPSI" = expectedPSI, "mutatedPSI" = mutatedPSI, row.names = NULL)
					}
					tempDf <- c("factor" = factor, "cancer" = cancer, "regressionCoefficient" = regressionCoefficient, "fdr" = fdr, "maxDistValue" = maxDistvalue, 
								minExpectation = minExp, maxExpectation = maxExp, "expectedPSI" = expectedPSI, "deviationPSI" = deviationPSI, "mutatedPSI" = mutatedPSI)
						
					mutatedPSIdf <- rbind(mutatedPSIdf, tempDf)
				}
				
				
			}
		}
		
	}
}


if (exonWise) {
	mutatedPSIdf <- data.frame(mutatedPSIdf)
	mutatedPSIdf$medianExpectedPSI <- as.numeric(as.character(mutatedPSIdf$medianExpectedPSI))
	mutatedPSIdf$medianMutatedPSI <- as.numeric(as.character(mutatedPSIdf$medianMutatedPSI))
	mutatedPSIdf$deltaOfMedians <- as.numeric(as.character(mutatedPSIdf$deltaOfMedians))
	mutatedPSIdf$regressionCoefficient <- as.numeric(as.character(mutatedPSIdf$regressionCoefficient))
	mutatedPSIdf$fractionIncreased <- as.numeric(as.character(mutatedPSIdf$fractionIncreased))
	mutatedPSIdf$fractionDecreased <- as.numeric(as.character(mutatedPSIdf$fractionDecreased))
	mutatedPSIdf$foldChange <- as.numeric(as.character(mutatedPSIdf$foldChange))
	#mutatedPSIdf$foldIncrease <- log2(mutatedPSIdf$fractionIncreased / mutatedPSIdf$fractionDecreased)
	
	ggplot(data = mutatedPSIdf, aes(x = foldChange, y = regressionCoefficient)) + geom_point(aes(color = factor(cancer))) + geom_smooth(method = "lm") + stat_cor()
	
	
	forplotDf <- mutatedPSIdf
	effectSize <- 0.05
	forplotDf$effectClass <- "NoChange"
	indexes <- forplotDf$deltaPSI < -effectSize
	forplotDf$effectClass[indexes] <- "Decreased"
	indexes <- forplotDf$deltaPSI > effectSize
	forplotDf$effectClass[indexes] <- "Increased"
	forplotDf <- forplotDf[forplotDf$effectClass != "NoChange", ]
	
	
	forplotDf <- mutatedPSIdf
	effectSize <- 0.58
	forplotDf$effectClass <- "NoChange"
	indexes <- forplotDf$foldChange < -effectSize
	forplotDf$effectClass[indexes] <- "Decreased"
	indexes <- forplotDf$foldChange > effectSize
	forplotDf$effectClass[indexes] <- "Increased"
	forplotDf <- forplotDf[forplotDf$effectClass != "NoChange", ]
	
	#forplotDf$effectClass <- ifelse(forplotDf$fractionDecreased > 0.1, "Decreased", "Increased")


	ggplot(data = forplotDf, aes(x = effectClass, y = regressionCoefficient)) + geom_boxplot() + facet_wrap(~cancer, scale = "free")
	table(forplotDf$effectClass, forplotDf$cancer)
	
}

if (!exonWise) {
	mutatedPSIdf <- data.frame(mutatedPSIdf)
	mutatedPSIdf$mutatedPSI <- as.numeric(as.character(mutatedPSIdf$mutatedPSI))
	mutatedPSIdf$expectedPSI <- as.numeric(as.character(mutatedPSIdf$expectedPSI))
	mutatedPSIdf$deviationPSI <- as.numeric(as.character(mutatedPSIdf$deviationPSI))


	mutatedPSIdf$deltaPSI <- mutatedPSIdf$mutatedPSI - mutatedPSIdf$expectedPSI
	mutatedPSIdf$shift <- mutatedPSIdf$deltaPSI / mutatedPSIdf$deviationPSI


	#save(mutatedPSIdf, file = "mutationAnalysisOfSplicingFactorsUsingSpecificEvents.Rda")
	#save(mutatedPSIdf, file = "mutationAnalysisOfSplicingFactorsUsingCommonEvents.Rda")

	plot(density(mutatedPSIdf$deviationPSI, na.rm = T))

	sdFilter <- 0.1
	forplotDf <- mutatedPSIdf[mutatedPSIdf$deviationPSI < sdFilter, ] %>% na.omit()

	forplotDf$regressionCoefficient <- as.numeric(as.character(forplotDf$regressionCoefficient))
	ggplot(data = forplotDf, aes(x = deltaPSI, y = regressionCoefficient)) + geom_point(aes(color = factor(cancer))) + geom_smooth(method = "lm") + stat_cor()



	#forplotDf <- mutatedPSIdf  %>% na.omit()
	#forplotDf$regressionCoefficient <- as.numeric(as.character(forplotDf$regressionCoefficient))
	#ggplot(data = forplotDf, aes(x = deltaPSI, y = regressionCoefficient)) + geom_point(aes(color = factor(cancer))) + geom_smooth(method = "lm") + stat_cor()

	effectSize <- 0.05
	forplotDf$effectClass <- "NoChange"
	indexes <- forplotDf$deltaPSI < -effectSize
	forplotDf$effectClass[indexes] <- "Decreased"
	indexes <- forplotDf$deltaPSI > effectSize
	forplotDf$effectClass[indexes] <- "Increased"
	forplotDf <- forplotDf[forplotDf$effectClass != "NoChange", ]


	#ggplot(data = forplotDf, aes(x = effectClass, y = regressionCoefficient)) + geom_boxplot() + facet_wrap(~cancer, scale = "free")
	plotFileName <- "mutationAnalysisOfSplicingFactorsUsingSpecificEvents.svg"
	ggplot(data = forplotDf, aes(x = effectClass, y = regressionCoefficient, fill = effectClass)) + geom_boxplot() + theme_bw() + ggTheme + 
		scale_fill_brewer(palette="RdBu") + stat_compare_means(method = "wilcox") + xlab("Embryonic splicing upon factor deletion")
	ggsave(file = plotFileName, height = 4, width = 4.5)
	table(forplotDf$effectClass, forplotDf$cancer)
	
	
	#ggplot(data = forplotDf, aes(x = deltaPSI, y = regressionCoefficient)) + geom_point(aes(color = factor(cancer))) + geom_smooth(method = "lm") + stat_cor()
	#ggplot(data = forplotDf, aes(x = deltaPSI, y = regressionCoefficient)) + geom_point() + geom_smooth() + stat_cor()
	
}



ggTheme <- theme(axis.text.x = element_text(size = 14, angle = 15, vjust =1, hjust = 1), legend.position = "none", #face = "bold"),
	axis.text.y = element_text(size = 14, face = "bold"), legend.title = element_blank(), plot.title = element_text(size = 10),
	axis.title = element_text(size = 14, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 1))


splicingFactors <- allFactors
allFactors <- tfList %>% unlist() %>% unique()

#############

mutatedPSIdf <- NULL
analysis <- "SFwise"
useSpecific <- T;

#load("commonEventsAndFactorsInLiverBrainKidney.Rda")

for (cancer in cancerType) {
	#mutationDf <- mutationDataList[[cancer]] %>% .[.$Variant_Classification == "Nonsense_Mutation" | .$SIFT == "deleterious", ]
	#mutationDf <- mutationDataList[[cancer]] %>% .[.$Variant_Classification == "Nonsense_Mutation", ]
	mutationDf1 <- mutationDataList[[cancer]] %>% .[.$Variant_Classification %in% mutationClasses, ]
	mutationDf2 <- mutationDataList[[cancer]] %>% .[grep(mutationConsequence, .$Consequence), ]
	
	mutationDf <- rbind(mutationDf1, mutationDf2) %>% unique()
	
	mutationDf$Tumor_Sample_Barcode <- gsub("\\-\\w+\\-\\w+\\-\\w+$", "", mutationDf$Tumor_Sample_Barcode)
	mutationDf$Tumor_Sample_Barcode <- gsub("[A-Z]$", "", mutationDf$Tumor_Sample_Barcode)
	mutationDf <- reshape2::dcast(mutationDf, Hugo_Symbol ~ Tumor_Sample_Barcode, fun.aggregate = length)
	
	mutationDf <- mutationDf %>% subset(., Hugo_Symbol %in% allFactors)
	mutatedSamples <- mutationDf[ ,-1] %>% colSums %>% .[. > 0] %>% names()
	
	cancerGeneTPM <- cancerExpressionList[[cancer]]
	
	rownames(cancerGeneTPM) <- gsub("-", "\\.", rownames(cancerGeneTPM))
	intactSamples <- colnames(cancerGeneTPM)[!(colnames(cancerGeneTPM) %in% mutatedSamples)]
	
	if (useSpecific) {
		events <- readoutList[[cancer]]
		readout <- cancerPSIlist[[cancer]] %>% .[rownames(.) %in% events, ] %>% apply(., 2, median, na.rm = T, na.action = na.pass)
	}else {
		readout <- cancerPSIlist[[cancer]] %>% .[rownames(.) %in% commonEvents, ] %>% apply(., 2, median, na.rm = T, na.action = na.pass)
	}
	
	
	for (index in 1:nrow(mutationDf)) {
		factor = mutationDf[index, 'Hugo_Symbol']
		#factorType = mutationDf[index, 'factorType']
		tempDf <- mutationDf[ ,-(1:2)] %>% data.matrix
		mutations <- tempDf[index, ] %>% unlist
	
		#intactSamples <- mutations[mutations == 0] %>% names
		mutatedSamples <- mutations[mutations > 0] %>% names
	
		if (analysis == "TFs") {
			for (mutatedSample in mutatedSamples) {
			
				presentInboth <- sum(names(cancerGeneTPM) %in% mutatedSample)
				if(presentInboth > 0) {
					desiredExpression <- cancerGeneTPM[rownames(cancerGeneTPM) %in% factor ,names(cancerGeneTPM) %in% mutatedSample] #%>% median(., na.rm = T, na.action = na.pass)
			
					cutoff <- (desiredExpression * 10) / 100
					allowedrange <- c(desiredExpression - cutoff, desiredExpression + cutoff)
			
					intactExpression <- cancerGeneTPM[rownames(cancerGeneTPM) %in% factor ,names(cancerGeneTPM) %in% intactSamples] %>% unlist()
			
					matchedExpressionSamples <- between(intactExpression, allowedrange[1], allowedrange[2])
					matchedExpressionSamples <- intactExpression[matchedExpressionSamples] %>% names()
					matchedExpression <- intactExpression[names(intactExpression) %in% matchedExpressionSamples] %>% median
					expectedPSI <- posPSI[names(posPSI) %in% matchedExpressionSamples] %>% median
					mutatedPSI <- posPSI[names(posPSI) %in% mutatedSample] %>% median
				
					factorExpression <- criticalExpression[names(criticalExpression) %in% matchedExpressionSamples] %>% median
					mutatedFactorExpression <- criticalExpression[names(criticalExpression) %in% mutatedSample] %>% median
				
					tempDf <- c("cancer" = unname(cancer), "TF" = factor, "TFexpression" = desiredExpression, "matchedTFexpression" = matchedExpression, 
							"splicingFactorIntact"= factorExpression, "splicingFactorMutated" = mutatedFactorExpression, "expectedPSI" = expectedPSI, "mutatedPSI" = mutatedPSI)
					mutatedPSIdf <- rbind(mutatedPSIdf, tempDf)
				}
			
			}
		}
		
		if (analysis == "SFwise") {
			for (mutatedSample in mutatedSamples) {
			
				presentInboth <- sum(names(cancerGeneTPM) %in% mutatedSample)
				if(presentInboth > 0) {
					desiredExpression <- cancerGeneTPM[rownames(cancerGeneTPM) %in% factor ,names(cancerGeneTPM) %in% mutatedSample] #%>% median(., na.rm = T, na.action = na.pass)
			
					cutoff <- (desiredExpression * 10) / 100
					allowedrange <- c(desiredExpression - cutoff, desiredExpression + cutoff)
			
					intactExpression <- cancerGeneTPM[rownames(cancerGeneTPM) %in% factor ,names(cancerGeneTPM) %in% intactSamples] %>% unlist()
			
					matchedExpressionSamples <- between(intactExpression, allowedrange[1], allowedrange[2])
					matchedExpressionSamples <- intactExpression[matchedExpressionSamples] %>% names()
					matchedExpression <- intactExpression[names(intactExpression) %in% matchedExpressionSamples] %>% median
					
					expectedPSI <- readout[names(readout) %in% matchedExpressionSamples] %>% median
					mutatedPSI <- readout[names(readout) %in% mutatedSample] %>% median
				
					splicingFactorExpression <- cancerGeneTPM[rownames(cancerGeneTPM) %in% splicingFactors ,names(cancerGeneTPM) %in% matchedExpressionSamples]
					medianExpression <- apply(splicingFactorExpression, 1, median, na.rm = T)
					deviation <- log(splicingFactorExpression + 1) %>% apply(., 1, sd, na.rm = T)
					
					mutatedSFexpression <- cancerGeneTPM[rownames(cancerGeneTPM) %in% splicingFactors ,names(cancerGeneTPM) %in% mutatedSample]
					
					foldChange <- log2(mutatedSFexpression / medianExpression)
					
					#tempDf <- merge(foldChange, coeffcientsDf, by = 0) %>% .[, c("Row.names", "x", cancer)] %>% set_names(c("factor", "coefficient", "cancer"))
					
					tempDf <- data.frame("cancer" = unname(cancer), "TF" = factor, "splicingFactor" = names(medianExpression), "TFexpression" = desiredExpression, "matchedTFexpression" = matchedExpression, 
							"splicingFactorIntact"= medianExpression, "splicingFactorDeviation" = deviation, "splicingFactorMutated" = mutatedSFexpression, "expectedPSI" = expectedPSI, "mutatedPSI" = mutatedPSI)
							
					mutatedPSIdf <- rbind(mutatedPSIdf, tempDf)
				}
			}
		}
	}
}
	




####

if (analysis == "TF") {
	mutatedPSIdf <- data.frame(mutatedPSIdf)
	mutatedPSIdf$mutatedPSI <- as.numeric(as.character(mutatedPSIdf$mutatedPSI))
	mutatedPSIdf$expectedPSI <- as.numeric(as.character(mutatedPSIdf$expectedPSI))
	mutatedPSIdf$splicingFactorMutated <- as.numeric(as.character(mutatedPSIdf$splicingFactorMutated))
	mutatedPSIdf$splicingFactorIntact <- as.numeric(as.character(mutatedPSIdf$splicingFactorIntact))


	mutatedPSIdf$deltaPSI <- (mutatedPSIdf$mutatedPSI) - mutatedPSIdf$expectedPSI
	mutatedPSIdf$splicingFactorFC <- log2((mutatedPSIdf$splicingFactorMutated / mutatedPSIdf$splicingFactorIntact))

	ggplot(data = mutatedPSIdf, aes(x = TF, y = splicingFactorFC)) + geom_boxplot()
	ggplot(data = mutatedPSIdf, aes(x = TF, y = deltaPSI)) + geom_boxplot()

	
}


if (analysis == "SFwise") {
	mutatedPSIdf <- data.frame(mutatedPSIdf)
	mutatedPSIdf$mutatedPSI <- as.numeric(as.character(mutatedPSIdf$mutatedPSI))
	mutatedPSIdf$expectedPSI <- as.numeric(as.character(mutatedPSIdf$expectedPSI))
	mutatedPSIdf$splicingFactorMutated <- as.numeric(as.character(mutatedPSIdf$splicingFactorMutated))
	mutatedPSIdf$splicingFactorIntact <- as.numeric(as.character(mutatedPSIdf$splicingFactorIntact))


	mutatedPSIdf$deltaPSI <- (mutatedPSIdf$mutatedPSI) - mutatedPSIdf$expectedPSI
	mutatedPSIdf$splicingFactorFC <- log2((mutatedPSIdf$splicingFactorMutated / mutatedPSIdf$splicingFactorIntact))
	
	tempDf <- coeffcientsDf %>% mutate(., splicingFactor = rownames(.)) %>% reshape2::melt() %>% set_colnames(c("splicingFactor", "cancer", "coefficient"))
	
	mutatedPSIdf <- merge(mutatedPSIdf, tempDf, by.x = c("cancer", "splicingFactor"),  by.y = c("cancer", "splicingFactor"))
	
	dataForPlot <- mutatedPSIdf[mutatedPSIdf$TF != "TP53", ]
	ggplot(data = dataForPlot, aes(x = deltaPSI, y = splicingFactorFC)) + geom_point() + facet_wrap(~cancer)
	#ggplot(data = mutatedPSIdf, aes(x = TF, y = splicingFactorFC)) + geom_boxplot()
	#ggplot(data = mutatedPSIdf, aes(x = TF, y = deltaPSI)) + geom_boxplot()
}



#############use this if euclidian distances should be weithed###########
####from https://www.statology.org/euclidean-distance-in-r/######
weightedDistance <- function(A, B, W) {
	#wSign <- sign(W)
	coefSums <- sum(W * (A - B)^2)
	coefSign <- sign(coefSums)
	coefSign * sqrt(abs(coefSums))
}

#W = combinedLoadings$coefficients
#W = combinedLoadings$coefficients * ifelse(combinedLoadings$FDR < 0.05, 1, 0)

#W = (1 /combinedLoadings$coefficients) %>% set_names(rownames(combinedLoadings))
W = (1 / combinedLoadings$coefficients) * ifelse(combinedLoadings$FDR < 0.05, 1, 0) %>% set_names(rownames(combinedLoadings))
W <- rescale(W, to = c(0, 100))
desiredFactors <- rownames(combinedLoadings)[rownames(combinedLoadings) %in% names(W)]
#W = (0.01 /combinedLoadings$coefficients) * ifelse(combinedLoadings$FDR < 0.05, 1, 0.1)


#predictors <- cancerGeneTPM %>% .[rownames(.) %in% criticalFactors, ]
#predictors <- cancerGeneTPM %>% .[rownames(.) %in% rownames(combinedLoadings), ]
predictors <- cancerGeneTPM %>% .[rownames(.) %in% desiredFactors, ]
predictors <- log2(predictors + 1)

weightDistDf <- list();
outputIndex = 0

for (index1 in 1:ncol(predictors)) {
	A <- predictors[, index1]
	for (index2 in 1:ncol(predictors)) {
		outputIndex <- outputIndex + 1
		B <- predictors[, index2]
		weightDist <- weightedDistance(A, B, W)
		temp <- data.frame(Var1 = colnames(predictors)[index1], Var2 = colnames(predictors)[index2], value = weightDist)
		weightDistDf[[outputIndex]] <- temp
		print(outputIndex)
	}
}
weightDistDf <- do.call("rbind", weightDistDf)





###################use this for overall mutation load vs splicing######################
mutationLoad <- mutationDf[ ,-1] %>% group_by(factorType) %>% summarise(across(everything(), sum), .groups = 'drop')  %>% 
			data.frame(check.names = F) %>% set_rownames(.$factorType) %>% .[ ,-1]  %>% t() %>% data.frame(check.names = F)

names(mutationLoad) <- gsub("\\w+\\.", "", names(mutationLoad))
names(posPSI) <- gsub("\\w+\\.", "", names(posPSI))
names(negPSI) <- gsub("\\w+\\.", "", names(negPSI))
names(criticalExpression) <- gsub("\\w+\\.", "", names(criticalExpression))
names(otherExpression) <- gsub("\\w+\\.", "", names(otherExpression))


mutationLoadDf <- merge(mutationLoad, posPSI, by = 0) %>% data.frame() %>% set_names(c("sample", "Critical", "NotCritical", "inclusionLevel"))
#mutationLoadDf <- merge(mutationLoad, negPSI, by = 0) %>% data.frame() %>% set_names(c("sample", "Critical", "Not Critical", "inclusionLevel"))

#mutationLoadDf <- merge(mutationLoadDf, criticalExpression, by.x = "sample", by.y = 0) %>% set_colnames(c("sample", "Critical", "Not Critical", "inclusionLevel", "criticalExpression"))
#mutationLoadDf <- merge(mutationLoadDf, otherExpression, by.x = "sample", by.y = 0) %>% set_colnames(c("sample", "Critical", "NotCritical", "inclusionLevel", "criticalExpression", "NotCriticalExpression"))
mutationLoadDf <- reshape2::melt(mutationLoadDf, id.vars = c("sample", "inclusionLevel")) 

colnames(mutationLoadDf) <- c("sample", "inclusionLevel", "factorType", "mutationLoad")
mutationLoadDf$sampleType <- ifelse(mutationLoadDf$mutationLoad == 0, "Not Mutated", "Mutated") %>% factor
mutationLoadDf$loadFactors <- factor(mutationLoadDf$mutationLoad)
#my_comparisons <- list(c("0", "1"), c("0", "11"))

ggplot(data = mutationLoadDf, aes(x = loadFactors, y = inclusionLevel)) + geom_boxplot() + facet_wrap(~factorType)

#subData <- mutationLoadDf[mutationLoadDf$criticalExpression > median(mutationLoadDf$criticalExpression), ]
#ggplot(data = subData, aes(x = loadFactors, y = inclusionLevel)) + geom_boxplot()





######maybe obsolete#############
#posPSI <- cancerPSI[rownames(cancerPSI) %in% embPos, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
#negPSI <- cancerPSI[rownames(cancerPSI) %in% embNeg, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)



ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	legend.title = element_blank(), legend.background = element_blank(),
	axis.title = element_text(size = 12, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	


################Analysis of mutation load vs embryonic splicing##############
rownames(combinedLoadings) <- gsub("\\.", "-", rownames(combinedLoadings))
criticalFactors <- combinedLoadings[combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0.0001, ] %>% rownames() #%>% GeneID2entrez(.)
OtherFactors <- rownames(combinedLoadings) #%>% GeneID2entrez(.)
OtherFactors <- OtherFactors[!(OtherFactors %in% criticalFactors)]
####select again top 100 genes####
criticalFactors <- combinedLoadings %>% .[.$FDR <= 0.05, ] %>% .[order(.$coefficients, decreasing = T), ] %>% rownames %>% .[1:100] #%>% GeneID2entrez(.)

criticalExpression <- cancerGeneTPM %>% .[rownames(.) %in% criticalFactors, ] %>% apply(., 2, median, )
otherExpression <- cancerGeneTPM %>% .[rownames(.) %in% OtherFactors, ] %>% apply(., 2, median, )




#mutationDf <- mutationDf %>% subset(., Hugo_Symbol %in% c(criticalFactors, OtherFactors))
#mutationDf <- mutationDf %>% mutate(., factorType = ifelse(.$Hugo_Symbol %in% criticalFactors, "Critical", "NotCritical"), .before = 2)


load(glue("TFbindingEnrichmentIn1kbPromotersUsingTFEA.ChIPremapDataIn{devTissue}SplicingFactors.Rda"))
topTFs <- rankedTFs[1:20, 'TF']

#mutationDf <- mutationSubData %>% .[.$Variant_Classification == "Nonsense_Mutation" | .$SIFT == "deleterious", ]
#mutationDf <- mutationSubData %>% .[.$Variant_Classification == "Nonsense_Mutation", ]
mutationDf <- mutationSubData %>% .[.$Variant_Classification == "Nonsense_Mutation" | .$SIFT == "deleterious" | .$PolyPhen == "probably_damaging", ]
#mutationDf <- mutationSubData %>% .[.$Variant_Classification == "Nonsense_Mutation" | (.$SIFT == "deleterious" & .$PolyPhen == "probably_damaging"), ]
mutationDf$Tumor_Sample_Barcode <- gsub("\\-\\w+\\-\\w+\\-\\w+$", "", mutationDf$Tumor_Sample_Barcode)
mutationDf$Tumor_Sample_Barcode <- gsub("[A-Z]$", "", mutationDf$Tumor_Sample_Barcode)

mutationDf <- reshape2::dcast(mutationDf, Hugo_Symbol ~ Tumor_Sample_Barcode)
mutationDf <- mutationDf %>% subset(., Hugo_Symbol %in% c(topTFs))
#mutationDf <- mutationDf %>% mutate(., factorType = ifelse(.$Hugo_Symbol %in% criticalFactors, "Critical", "NotCritical"), .before = 2)


intactSamples <- mutationDf[ ,-(1:2)] %>% colSums %>% .[. == 0] %>% names()
#intactExpression <- criticalExpression[names(criticalExpression) %in% intactSamples]


