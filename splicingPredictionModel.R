#library(clusterProfiler)
library(dplyr); library(magrittr); library(ensembldb); library(EnsDb.Hsapiens.v86)
library(gtools); library(preprocessCore); library(reshape2); library(ggplot2)
library(clusterProfiler); library(msigdbr); library(gplots); library(FactoMineR); library(factoextra)
library(glue); library(pls); library(TCGAbiolinks)
orgDb <- "org.Hs.eg.db"; EnsdB <- EnsDb.Hsapiens.v86

Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),

GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))

processEpxressionFile <- function(dataFrame, Transcripts, replicates) {
	rownames(dataFrame) <- gsub("\\.\\d+", "", rownames(dataFrame))
	matches <- match(rownames(dataFrame), Transcripts$tx_id)
	Genes <- Transcripts[matches, 'gene_name']
	dataFrame <- cbind("GeneName" = Genes, dataFrame)
	GeneLevelTPM <- aggregate(.~GeneName, data = dataFrame, sum, na.rm = T, na.action = na.pass) 
	GeneLevelTPM <- GeneLevelTPM[ ,-1] %>% as.matrix() %>% normalize.quantiles() %>% 
			set_rownames(GeneLevelTPM$GeneName) %>% 
			set_colnames(colnames(GeneLevelTPM)[-1]) %>% data.frame(check.names = F)
		
	if (replicates) {
		GeneLevelTPM <- GeneLevelTPM %>% t() %>% data.frame() %>%
				mutate(Stage = gsub("_Rep\\d+", "", colnames(GeneLevelTPM))) %>% 
				aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
				set_rownames(.$Stage) %>% .[,-1] %>% t() %>% data.frame(check.names = F) %>%
				set_rownames(rownames(GeneLevelTPM))
	}	
	
	#GeneLevelTPM <- GeneLevelTPM[ ,mixedsort(names(GeneLevelTPM))]
	return(GeneLevelTPM)
}

processPSIfile <- function(psiFile) {
	data.frame(t(psiFile)) %>% 
	mutate(Stage = gsub("_Rep\\d+", "", names(psiFile))) %>% 
	aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
	set_rownames(.$Stage) %>% 
	.[,-1] %>% t() %>% data.frame(check.names = F) %>%
	set_rownames(rownames(psiFile))
}

exonGeneNames <- function(exonList, GenesDf) {
	exonGenes <-  exonList %>% as.character() %>% gsub("\\..+", "", .)
	indexes <- match(exonGenes, GenesDf[ ,'GeneID'])
	exonGenes <- GenesDf[indexes, c('GeneName')] %>% as.character()
	geneType <- GenesDf[indexes, c('GeneType')] %>% as.character()
	return(cbind(exonGenes, geneType))
}

tcgaStaging <- function(cancerType, samples) {
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
	
	samples <- gsub("-\\w+$", "", samples)
	clinical <- clinical[match(samples, clinical$Barcode), ]

	early = c("Stage I", "Stage II"); late = c("Stage III", "Stage IV")
	earlySamples <- clinical[clinical$Stage %in% early, 'Barcode'] %>% as.character()
	lateSamples <- clinical[clinical$Stage %in% late, 'Barcode'] %>% as.character()
	stageDf <- rbind(data.frame(sampleName = earlySamples, sampleType = "Early Cancer"), data.frame(sampleName = lateSamples, sampleType = "Late Cancer"))
	
	return(stageDf)
}

loadingFunction <- function(model) {
	load <- with(model, unclass(loadings))
	aload <- abs(load) ## save absolute values
	tmpDf <- sweep(aload, 2, colSums(aload), "/")
	return(tmpDf)
}

medianImpute <- function(x){
	ifelse(is.na(x), median(x, na.rm = T), x)
}

CorrelationFunction <- function(index, PSIfile, pathwayFile) {
	exonUsage <- PSIfile[, index] %>% unlist()
	corList <- list()
	pvalList <- list()
	outputList <- list()
	for(pathwayIndex in 1:ncol(pathwayFile)) {
		pathwayUsage <- pathwayFile[ ,pathwayIndex]
		correlations <- cor.test(exonUsage, pathwayUsage)
		corList[pathwayIndex] <- correlations$estimate
		pvalList[pathwayIndex] <-  correlations$p.value
	}
	#outputList[['corr']] <- unlist(corList)
	#outputList[['pval']] <- unlist(pvalList)
	#return(outputList)
	cbind(correlation = unlist(corList), pvalue = unlist(pvalList))
}

binTestingFunction <- function(PSIfile, nbins) {
	tempDf <- apply(PSIfile, 1, function(x) {y = ((x*100)/nbins) %>% floor(); ifelse(y <= 4, y+1, y-1)})
	lengths <- apply(tempDf, 2, function(x) {length(unique(x))})
	indexes1 <- lengths > 1
	indexes2 <- tempDf[ ,indexes1] %>% apply(., 2, function(x) {ifelse(length(unique(x)) == 2, {x1 = unique(x)[1]; x2 = unique(x)[2]; ifelse(abs(x2 - x1) == 1, F, T)}, T)})		
	tempDf[ ,indexes1] %>% .[, indexes2] %>% t() %>% rownames()
}



#pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/Normal/Kallisto/"
#tissueCancerFile <- read.table("../Cancer/Kallisto/cancerAndTissuePairs", sep = "\t", header = T)

splicingFactors <- read.table("/Users/singha30/DataAndAnalysis/ExonicRepeatsInCancer/SplicingFactorsGO-0008380fromAmigo", sep = "\t", header = F)

tissueCancerFile <- read.table("../Cancer/Kallisto/cancerAndTissuePairs", sep = "\t", header = T)
tissueTypes <- tissueCancerFile[ ,'devTissue'] %>% as.character() %>% unique()
pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/"


regressAll = F
fileIndex = 7
#for (fileIndex in 1:nrow(tissueCancerFile)) {
	tissueType <- tissueCancerFile[fileIndex ,'Tissue'] %>% as.character()
	cancerType <- tissueCancerFile[fileIndex ,'Cancer'] %>% as.character()
	devTissue <- tissueCancerFile[fileIndex ,'devTissue'] %>% as.character()
	
	load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))
	
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType != "Not defined", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	
	embPos <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Positive", 'Exon'] %>% as.character()
	embNeg <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Negative", 'Exon'] %>% as.character()
		
	devFile <- glue("{pathToExpressionFile}Normal/Kallisto/development/PSIvalues/SuppaPSIvaluesUsingKallistoIn{devTissue}.psi")
	devPSIfile <- read.table(devFile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)
	
	NAsums <- apply(devPSIfile,1, function(x) sum(is.na(x)))
	NAthreshold <- ceiling(ncol(devPSIfile) / 2)
	indexes <- which(NAsums < NAthreshold)
	devPSIfile <- devPSIfile[indexes, ]

	devPSIfile <- processPSIfile(psiFile = devPSIfile)
	devPSIfile <- devPSIfile[, mixedsort(names(devPSIfile))]
	
	posPSI <- devPSIfile[rownames(devPSIfile) %in% embPos, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
	negPSI <- devPSIfile[rownames(devPSIfile) %in% embNeg, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
	
	#variableMatrix <- data.frame(posPSI = posPSI, negPSI = negPSI) %>% as.matrix()
	#variableMatrix <- data.frame(negPSI = negPSI) %>% as.matrix()
	variableMatrix <- data.frame(posPSI = posPSI) %>% as.matrix()
	#variableMatrix <- scale(variableMatrix, center = T, scale = F)
	
	if (regressAll) {
		exonsToUse <- binTestingFunction(PSIfile = devPSIfile, nbins = 20)
		variableMatrix <- devPSIfile[rownames(devPSIfile) %in% exonsToUse, ]
		variableMatrix <- apply(variableMatrix, 1, function(x) medianImpute(x))
	}
	
	devFile <- glue("{pathToExpressionFile}Normal/Kallisto/development/Expression/KallistoTPMvaluesForTranscriptExpressionIn{devTissue}")
	cancerExpression <- glue("{pathToExpressionFile}Cancer/Kallisto/TCGA/Expression/KallistoTPMvaluesForTranscriptExpressionIn{cancerType}")
	normalExpression <- glue("{pathToExpressionFile}Cancer/Kallisto/Gtex/Expression/KallistoTPMvaluesForTranscriptExpressionIn{tissueType}")
	
	devExpression <- read.table(paste(devFile), sep = "\t", header = T, check.names = F)
	cancerExpression <- read.table(cancerExpression, sep = "\t", header = T, check.names = F)
	normalExpression <- read.table(normalExpression, sep = "\t", header = T, check.names = F)
	
	names(cancerExpression) <- gsub("\\.", "-", names(cancerExpression))
	names(normalExpression) <- gsub("\\.", "-", names(normalExpression))
	
	cancerSamples <- TCGAquery_SampleTypes(names(cancerExpression), typesample = c("TP"))
	cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples]
	
	commonTranscripts <- rownames(cancerExpression)[rownames(cancerExpression) %in% rownames(normalExpression)] %>% .[. %in% rownames(devExpression)]
							
	cancerExpression <- cancerExpression[rownames(cancerExpression) %in% commonTranscripts, ] %>% .[order(rownames(.)), ]
	normalExpression <- normalExpression[rownames(normalExpression) %in% commonTranscripts, ] %>% .[order(rownames(.)), ]
	devExpression <- devExpression[rownames(devExpression) %in% commonTranscripts, ] %>% .[order(rownames(.)), ]
	
	combinedFile <- cbind(devExpression, cancerExpression, normalExpression)
	geneLevelTPM <- processEpxressionFile(dataFrame = combinedFile, Transcripts = Transcripts, replicates = T)

	cancerGeneTPM <- geneLevelTPM[ ,names(geneLevelTPM) %in% names(cancerExpression)]
	normalGeneTPM <- geneLevelTPM[ ,names(geneLevelTPM) %in% names(normalExpression)]
	devNames <- gsub("_Rep\\d+", "", colnames(devExpression))
	devGeneTPM <- geneLevelTPM[ ,names(geneLevelTPM) %in% devNames]
	devGeneTPM <- devGeneTPM[ ,mixedsort(names(devGeneTPM))]
	
	#devGeneTPM <- processEpxressionFile(dataFrame = devExpression, Transcripts = Transcripts, replicates = T)
	
	
	predictorMatrix <- devGeneTPM[rownames(devGeneTPM) %in% splicingFactors[,2], ] %>% t() %>% data.frame()
	testMatrix <- cancerGeneTPM[rownames(cancerGeneTPM) %in% splicingFactors[,2], ] %>% t() %>% data.frame()
	
	predictorMatrix <- log2(predictorMatrix + 1) %>% as.matrix()
	testMatrix <- log2(testMatrix + 1) %>% as.matrix()
	
	indexes <- apply(predictorMatrix, 2, sd, na.rm = T) > 0 & apply(testMatrix, 2, sd, na.rm = T) > 0
	predictorMatrix <- predictorMatrix[ ,indexes]
	testMatrix <- testMatrix[ ,indexes]
	
	set.seed(123)
	model <- plsr(variableMatrix ~ predictorMatrix, jackknife = T, validation = "LOO")
	
	
	
	####load the data to test the model####
	cancerPSIfile <- glue("{pathToExpressionFile}Cancer/Kallisto/TCGA/PSIvalues/SuppaPSIvaluesUsingKallistoIn{cancerType}.psi")
	cancerPSI <- read.table(cancerPSIfile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)
	names(cancerPSI) <- gsub("\\.", "-", names(cancerPSI))
	
	NAsums <- apply(cancerPSI,1, function(x) sum(is.na(x)))
	NAthreshold <- ceiling(ncol(cancerPSI) / 2)
	indexes <- which(NAsums < NAthreshold)
	cancerPSI <- cancerPSI[indexes, ]
	cancerPSI <- processPSIfile(psiFile = cancerPSI)
	cancerPSI <- cancerPSI[, names(cancerPSI) %in% cancerSamples]
	
	posPSI <- cancerPSI[rownames(cancerPSI) %in% embPos, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
	negPSI <- cancerPSI[rownames(cancerPSI) %in% embNeg, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
	
	#variableMatrix <- data.frame(posPSI = posPSI, negPSI = negPSI) %>% as.matrix()
	#testVariables <- data.frame(negPSI = negPSI) %>% as.matrix()
	testVariables <- data.frame(posPSI = posPSI) %>% as.matrix()
	#testVariables <- scale(testVariables, center = T, scale = F)
	predictedValues <- predict(model, ncomp = 2, newdata = testMatrix)
	plot(testVariables, predictedValues)
	
	contributionDf <- loadingFunction(model)
	coefficientsDf <- coef(model,1)
	scoresDf <- scores(model)
	
	R2object <- R2(model, ncomp = 2, estimate = "train")
	Yvar <- R2object$val
	
	jackTest <- jack.test(model, ncomp = 2, use.mean = TRUE)
	Pvalues <- jackTest$pvalues
	deviation <- jackTest$sd
	
	combinedLoadings <- data.frame(PC1.contribution = contributionDf[,1:2], coefficients = coefficientsDf[,1,1], pvalue = Pvalues[,1,1], standardDevation = deviation[,,1])
	combinedLoadings$FDR <- p.adjust(combinedLoadings$pvalue, method = "fdr")
	#plot(scoresDf[,1], variableMatrix[,1])
	save(model, combinedLoadings, file = glue("plsrModelForMedianPosEmbSplicingIn{devTissue}.rda"))
	#break
	#}




	######################
	
dataForplot <- data.frame(actualValue = testVariables[ ,1], predictedValue = predictedValues[,1,1])
ggplot(data = dataForplot, aes(x = actualValue, y = predictedValue)) + geom_point(col = "mediumpurple2") + theme_bw() + geom_smooth(method = "lm", se = F)