library(TCGAbiolinks); library(ensembldb); 
library(EnsDb.Hsapiens.v86); EnsdB <- EnsDb.Hsapiens.v86
library(magrittr); library(dplyr); library(TCGAbiolinks)
library(preprocessCore); library(glue)
####load DNA sequence if needed######
#dna <- ensembldb:::getGenomeTwoBitFile(EnsdB)


Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),

GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))

pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/Cancer/Kallisto/"
tissueCancerFile <- read.table("cancerAndTissuePairs", sep = "\t", header = T)


processEpxressionFile <- function(dataFrame, Transcripts) {
	rownames(dataFrame) <- gsub("\\.\\d+", "", rownames(dataFrame))
	matches <- match(rownames(dataFrame), Transcripts$tx_id)
	Genes <- Transcripts[matches, 'gene_name']
	dataFrame <- cbind("GeneName" = Genes, dataFrame)
	GeneLevelTPM <- aggregate(.~GeneName, data = dataFrame, sum, na.rm = T, na.action = na.pass) 
	GeneLevelTPM <- GeneLevelTPM[ ,-1] %>% as.matrix() %>% normalize.quantiles() %>% 
			set_rownames(GeneLevelTPM$GeneName) %>% 
			set_colnames(colnames(GeneLevelTPM)[-1]) %>% data.frame(check.names = F)
			
	return(GeneLevelTPM)
}

wilcoxTest <- function(index, condition1, condition2, allowedNA) {
	sample1 <- condition1[index, ] %>% unlist()
	sample2 <- condition2[index, ] %>% unlist()
	pvalue = NA; effectSize = NA
	normalMedian = NA; cancerMedian = NA
	exon <- rownames(condition1)[index]
	normalNAcutoff <- (length(sample1) * (allowedNA / 100)) %>% floor()
	cancerNAcutoff <- (length(sample2) * (allowedNA / 100)) %>% floor()
	if((sum(is.na(sample1)) <= normalNAcutoff) & (sum(is.na(sample2)) <= cancerNAcutoff)) {
		res <- wilcox.test(sample1, sample2)
		pvalue <- res$p.value
		normalMedian <- median(sample1, na.rm = T, na.action = na.pass)
		cancerMedian <- median(sample2, na.rm = T, na.action = na.pass)
		effectSize <- cancerMedian - normalMedian
	}
	cbind(exon, normalMedian, cancerMedian, effectSize, pvalue)
}

OutlierDetectionFunction <- function(index, condition1, condition2, allowedNA, numberOfSDs, normalCutOff, cancerCutOff) {
	sample1 <- condition1[index, ] %>% unlist()
	sample2 <- condition2[index, ] %>% unlist()
	
	normalMean = NA; cancerMean = NA; 
	normalDeviation = NA; eventType = NA
	
	exon <- rownames(condition1)[index]
	normalNAcutoff <- (length(sample1) * (allowedNA / 100)) %>% floor()
	cancerNAcutoff <- (length(sample2) * (allowedNA / 100)) %>% floor()
	if((sum(is.na(sample1)) <= normalNAcutoff) & (sum(is.na(sample2)) <= cancerNAcutoff)) {
		normalMean <- mean(sample1, na.rm = T, na.action = na.pass)
		cancerMean <- mean(sample2, na.rm = T, na.action = na.pass)
		normalDeviation <- sd(sample1, na.rm = T)
		upperLimit <- normalMean + (numberOfSDs * normalDeviation)
		lowerLimit <- normalMean - (numberOfSDs * normalDeviation)
		#cancerMedian <- median(sample2, na.rm = T, na.action = na.pass)
		#effectSize <- cancerMedian - normalMedian
		
		normalOutlierUpper <- sum(sample1 > upperLimit, na.rm = T)
		normalOutlierLower <- sum(sample1 < lowerLimit, na.rm = T)
		
		normalOutlierUpper <- (normalOutlierUpper / length(na.omit(sample1))) * 100
		normalOutlierLower <- (normalOutlierLower / length(na.omit(sample1))) * 100
		
		
		frequentIncrease <- sum(sample2 > upperLimit, na.rm = T)
		frequentDecrease <- sum(sample2 < lowerLimit, na.rm = T)
		
		percentIncrease <- (frequentIncrease / length(na.omit(sample2))) * 100
		percentDecrease <- (frequentDecrease / length(na.omit(sample2))) * 100
		
		
		if (normalOutlierUpper < normalCutOff) {
			if (percentIncrease >= cancerCutOff) {
				eventType = "Increase"
			}
		}
		
		if (normalOutlierLower < normalCutOff) {
			if (percentDecrease >= cancerCutOff) {
				eventType = "Decrease"
			}
		}
		
		if ((normalOutlierUpper < normalCutOff) & (normalOutlierLower < normalCutOff)) {
			if ((percentIncrease >= cancerCutOff) & (percentDecrease >= cancerCutOff)) {
				eventType = "TwoWays"
			}
		}
		
		if (is.na(eventType)) {
			eventType = "Other"
		}
	}
	cbind(exon, normalMean, cancerMean, normalDeviation, eventType)
}


wilcoxParseFunction <- function(file, nRows, nCols, colNames, effectSize, FDR, cancerStage) {
	unlist(file) %>% matrix(., nrow = nRows, ncol = nCols, byrow = T) %>% data.frame(.) %>%
	na.omit() %>% set_colnames(colNames) %>%
	mutate(medianDeltaPSI = as.numeric(as.character(.$medianDeltaPSI)), Pvalue = as.numeric(as.character(.$Pvalue))) %>%
	cbind(., exonGeneNames(exonList = .$Exon, GenesDf = GenesDf)) %>% set_colnames(c(colNames, "GeneName", "GeneType")) %>%
	DifferentialExons(DiffPSI = ., geneType = "protein_coding", effectSize = effectSize, FDR = FDR) %>% mutate(CancerStage = cancerStage)
}

outlierParseFunction <- function(file, nRows, nCols, colNames, cancerStage) {
	unlist(file) %>% matrix(., nrow = nRows, ncol = nCols, byrow = T) %>% data.frame(.) %>%
	na.omit() %>% set_colnames(colNames) %>%
	cbind(., exonGeneNames(exonList = .$Exon, GenesDf = GenesDf)) %>% set_colnames(c(colNames, "GeneName", "GeneType")) %>%
	.[.$eventType != "Other", ] %>% mutate(CancerStage = cancerStage)
}

exonGeneNames <- function(exonList, GenesDf) {
	exonGenes <-  exonList %>% as.character() %>% gsub("\\..+", "", .)
	indexes <- match(exonGenes, GenesDf[ ,'GeneID'])
	exonGenes <- GenesDf[indexes, c('GeneName')] %>% as.character()
	geneType <- GenesDf[indexes, c('GeneType')] %>% as.character()
	return(cbind(exonGenes, geneType))
}

DifferentialExons <- function(DiffPSI, geneType, effectSize, FDR) {
	if (geneType == "All") {
		DiffPSI <- DiffPSI
	} else {
		DiffPSI <- DiffPSI[DiffPSI$GeneType == geneType, ]
	}
	DiffPSI$FDR <- p.adjust(DiffPSI$Pvalue, method = "fdr")
	indexes <- abs(DiffPSI$medianDeltaPSI) >= effectSize & DiffPSI$FDR <= FDR
	DiffPSI[indexes, ]
}

paired = F; markerBasedStaging = F
DifferentialEventsDf <- NULL
sampleCountsDf <- NULL

test = "outlier"
for (fileIndex in 1:nrow(tissueCancerFile)) {
	tissueType <- tissueCancerFile[fileIndex ,'Tissue'] %>% as.character()
	cancerType <- tissueCancerFile[fileIndex ,'Cancer'] %>% as.character()
	
	#if(cancerType == "GBM" | cancerType == "OV") next #####becuase these cancer types didn't have cancer type information.
	normalFile <- glue("{pathToExpressionFile}Gtex/Expression/KallistoTPMvaluesForTranscriptExpressionIn{tissueType}")
	cancerFile <- glue("{pathToExpressionFile}TCGA/Expression/KallistoTPMvaluesForTranscriptExpressionIn{cancerType}")
	normalExpression <- read.table(paste(normalFile), sep = "\t", header = T, check.names = F)
	cancerExpression <- read.table(paste(cancerFile), sep = "\t", header = T, check.names = F)
	names(cancerExpression) <- gsub("\\.", "-", names(cancerExpression))
	names(normalExpression) <- gsub("\\.", "-", names(normalExpression))
	
	if (paired) {
		matchedSamples <- TCGAquery_MatchedCoupledSampleTypes(names(cancerExpression), c("NT", "TP"))
		offset = length(matchedSamples) / 2
		normalSamples <- matchedSamples[1:offset]
		cancerSamples <- matchedSamples[(offset + 1):(offset * 2)]
		normalExpression <- cancerExpression[, names(cancerExpression) %in% normalSamples] %>% .[ ,order(names(.))]
		cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples] %>% .[ ,order(names(.))]
	
	} else {
		cancerSamples <- TCGAquery_SampleTypes(names(cancerExpression), typesample = c("TP"))
		cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples]
	}
	
	normalGeneTPM <- processEpxressionFile(dataFrame = normalExpression, Transcripts = Transcripts) 
	cancerGeneTPM <- processEpxressionFile(dataFrame = cancerExpression, Transcripts = Transcripts)
	
	normalMedian <- apply(normalGeneTPM, 1, median, na.rm = T, na.action = na.pass)
	cancerMedian <- apply(cancerGeneTPM, 1, median, na.rm = T, na.action = na.pass)
	#######both these files should be having matched column names in same order##########
	
	normalFile <- glue("{pathToExpressionFile}Gtex/PSIvalues/SuppaPSIvaluesUsingKallistoIn{tissueType}.psi")
	cancerFile <- glue("{pathToExpressionFile}TCGA/PSIvalues/SuppaPSIvaluesUsingKallistoIn{cancerType}.psi")
	
	normalPSI <- read.table(paste(normalFile), sep = "\t", header = T, check.names = F, na.strings = c(NA, "nan"))
	cancerPSI <- read.table(paste(cancerFile), sep = "\t", header = T, check.names = F, na.strings = c(NA, "nan"))
	names(cancerPSI) <- gsub("\\.", "-", names(cancerPSI))
	names(normalPSI) <- gsub("\\.", "-", names(normalPSI))

	if (paired) {
		normalPSI <- cancerPSI[, match(normalSamples, names(cancerPSI))]  %>% .[ ,order(names(.))]
		cancerPSI <- cancerPSI[, match(cancerSamples, names(cancerPSI))]  %>% .[ ,order(names(.))]
		diffPSI <- cancerFile - normalFile
	}else {
		cancerSamples <- TCGAquery_SampleTypes(names(cancerPSI), typesample = c("TP"))
		cancerPSI <- cancerPSI[, names(cancerPSI) %in% cancerSamples]
	}
	
	
	commonExons <- rownames(normalPSI)[rownames(normalPSI) %in% rownames(cancerPSI)]
	normalPSI <- normalPSI[rownames(normalPSI) %in% commonExons, ]
	cancerPSI <- cancerPSI[rownames(cancerPSI) %in% commonExons, ]
	
	stageInformation <- "Not defined"
	clinical <- GDCquery_clinic(paste("TCGA-",cancerType, sep = ""), type = "clinical")
	if (is.null(clinical$ajcc_pathologic_stage)) {
		stageInformation = "Not defined"
		if (cancerType == "OV") {
			clinical <- data.frame(Barcode = clinical$submitter_id, Stage = clinical$figo_stage)#, MetaStasis = clinical$ajcc_pathologic_m)
			clinical$Stage <- gsub("[A-HJ-UW-Z]$", "", clinical$Stage)
			if (length(table(clinical$Stage)) > 0) {
				stageInformation = "defined"
			}
		}
	} else {
		clinical <- data.frame(Barcode = clinical$submitter_id, Stage = clinical$ajcc_pathologic_stage)#, MetaStasis = clinical$ajcc_pathologic_m)
		clinical$Stage <- gsub("[A-HJ-UW-Z]$", "", clinical$Stage)
		if (length(table(clinical$Stage)) > 0) {
			stageInformation = "defined"
		}
	}

	names(cancerPSI) <- gsub("-\\w+$", "", names(cancerPSI))
	samplesUsed <- names(cancerPSI)
	clinical <- clinical[match(samplesUsed, clinical$Barcode), ]

	early = c("Stage I", "Stage II"); late = c("Stage III", "Stage IV")
	earlySamples <- clinical[clinical$Stage %in% early, 'Barcode'] %>% as.character()
	lateSamples <- clinical[clinical$Stage %in% late, 'Barcode'] %>% as.character()
	
	earlyPSI <- cancerPSI[ ,names(cancerPSI) %in% earlySamples]
	latePSI <- cancerPSI[ ,names(cancerPSI) %in% lateSamples]
	
	#####clinical data#####
	if(markerBasedStaging) {
		clinical <- read.table(paste("hallmarksBasedStagingFor", cancerType, ".txt", sep = ""), sep = "\t", header = T)
		earlySamples <- clinical[clinical$stage == "Early", 'samples']
		lateSamples <- clinical[clinical$stage == "Late", 'samples']
	
		earlyPSI <- cancerPSI[ ,names(cancerPSI) %in% earlySamples]
		latePSI <- cancerPSI[ ,names(cancerPSI) %in% lateSamples]
	}
	
	sampleCounts <- c(cancer = cancerType, Early = ncol(earlyPSI), Late = ncol(latePSI))
	sampleCountsDf <- rbind(sampleCountsDf, sampleCounts)
	
	overallDiffPSI <- NULL; earlyDiffPSI <- NULL; lateDiffPSI <- NULL
	
	
	if (test == "outlier") {
		
		colNames <- c("Exon", "normalMean", "cancerMean", "normalDeviation", "eventType")
		
		overallDiffPSI <- sapply(1:nrow(cancerPSI), OutlierDetectionFunction, condition1 = normalPSI, condition2 = cancerPSI, allowedNA = 50, numberOfSDs = 2, normalCutOff = 5, cancerCutOff = 20) %>% 
				outlierParseFunction(file = ., nRows = nrow(cancerPSI), nCols = 5, colNames = colNames, cancerStage = "Overall")
				
		
		if (stageInformation == "defined")	{
			earlyDiffPSI <- sapply(1:nrow(earlyPSI), OutlierDetectionFunction, condition1 = normalPSI, condition2 = earlyPSI, allowedNA = 50, numberOfSDs = 2, normalCutOff = 5, cancerCutOff = 20) %>% 
					outlierParseFunction(file = ., nRows = nrow(earlyPSI), nCols = 5, colNames = colNames, cancerStage = "Early")
				
				
			lateDiffPSI <- sapply(1:nrow(latePSI), OutlierDetectionFunction, condition1 = normalPSI, condition2 = latePSI, allowedNA = 50, numberOfSDs = 2, normalCutOff = 5, cancerCutOff = 20) %>% 
					outlierParseFunction(file = ., nRows = nrow(latePSI), nCols = 5, colNames = colNames, cancerStage = "Late")	
		}
		
		FrequentEventsDf <- rbind(overallDiffPSI, earlyDiffPSI, lateDiffPSI) %>% mutate(CancerType = cancerType, TissueType = tissueType)
	}
	
	earlyIncreaseExons <- FrequentEventsDf$Exon[which(FrequentEventsDf$CancerStage == "Early" & FrequentEventsDf$eventType == "Increase")]
	earlyIncreaseIndexInLate <- which(FrequentEventsDf$Exon %in% earlyIncreaseExons & FrequentEventsDf$eventType == "Increase" & FrequentEventsDf$CancerStage == "Late")
	FrequentEventsDf$CancerStage[earlyIncreaseIndexInLate] <- "Early.Late"

	earlyDecreaseExons <- FrequentEventsDf$Exon[which(FrequentEventsDf$CancerStage == "Early" & FrequentEventsDf$eventType == "Decrease")]
	earlyDecreaseIndexInLate <- which(FrequentEventsDf$Exon %in% earlyDecreaseExons & FrequentEventsDf$eventType == "Decrease" & FrequentEventsDf$CancerStage == "Late")
	FrequentEventsDf$CancerStage[earlyDecreaseIndexInLate] <- "Early.Late"
	
	indexes <- match(FrequentEventsDf$GeneName, names(normalMedian))
	FrequentEventsDf$NormalExpression <- normalMedian[indexes]
	
	indexes <- match(FrequentEventsDf$GeneName, names(cancerMedian))
	FrequentEventsDf$CancerExpression <- cancerMedian[indexes]
	
	DifferentialEventsDf <- rbind(DifferentialEventsDf, FrequentEventsDf)
	#break
}

write.table(DifferentialEventsDf, glue("DifferentialEventsInCancerWithTCGAstagingUsing{test}TestFromKallistoTrancripts-Comprehensive-2SD.txt"), sep = "\t", quote = F, row.names = F)
write.table(sampleCountsDf, "SampleCountDataFrame-Comprehensive.txt", sep = "\t", quote = F, row.names = F)

