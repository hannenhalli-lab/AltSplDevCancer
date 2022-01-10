#library(clusterProfiler)
library(dplyr); library(magrittr); library(ensembldb); library(EnsDb.Hsapiens.v86)
library(preprocessCore); library(reshape2); library(ggplot2)
library(glue); library(pls); library(TCGAbiolinks)
orgDb <- "org.Hs.eg.db"; EnsdB <- EnsDb.Hsapiens.v86


source('helperFunctions.R')

Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),

GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))




getEvents <- function (AssociationClusters, ExonType, Association) {
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	ExonType <- paste(ExonType, Association, sep = ".")
	events <- AssociationClusters[AssociationClusters$ExonType == ExonType, 'Exon'] %>% as.character()
	return(events)
}



pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/"


devTissues <- c("Hindbrain", "Liver", "Kidney")
testCancers <- c(Hindbrain = "LGG", Liver = "LIHC", Kidney = "KIRP")

pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/"


cancerPredictionList <- list()
cancerPredictedList <- list()
for (cancerType in testCancers) {
	
	cancerExpression <- glue("{pathToExpressionFile}Cancer/Kallisto/TCGA/Expression/KallistoTPMvaluesForTranscriptExpressionIn{cancerType}")
	cancerExpression <- read.table(cancerExpression, sep = "\t", header = T, check.names = F)
	names(cancerExpression) <- gsub("\\.", "-", names(cancerExpression))
	cancerSamples <- TCGAquery_SampleTypes(names(cancerExpression), typesample = c("TP"))
	cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples]
	
	geneLevelTPM <- processEpxressionFile(dataFrame = cancerExpression, Transcripts = Transcripts, replicates = F)
	rownames(geneLevelTPM) <- gsub("-", "\\.", rownames(geneLevelTPM))
	
	cancerPSIfile <- glue("{pathToExpressionFile}Cancer/Kallisto/TCGA/PSIvalues/SuppaPSIvaluesUsingKallistoIn{cancerType}.psi")
	cancerPSI <- read.table(cancerPSIfile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)
	names(cancerPSI) <- gsub("\\.", "-", names(cancerPSI))
	cancerPSI <- cancerPSI[, names(cancerPSI) %in% cancerSamples]
	
	NAsums <- apply(cancerPSI,1, function(x) sum(is.na(x)))
	NAthreshold <- ceiling(ncol(cancerPSI) / 2)
	indexes <- which(NAsums < NAthreshold)
	cancerPSI <- cancerPSI[indexes, ]
	#cancerPSI <- processPSIfile(psiFile = cancerPSI) #### incluse this if replicates
	
	
	devTissue <- testCancers[testCancers %in% cancerType]  %>% names
	load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))
	embPos <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Positive")
	posPSI <- cancerPSI[rownames(cancerPSI) %in% embPos, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
	
	
	tissueList1 <- list(); tissueList2 <- list()
	for (devTissue in devTissues) {
		
		load(glue("plsrModelForMedianPosEmbSplicingIn{devTissue}.rda"))
		testMatrix <- geneLevelTPM[rownames(geneLevelTPM) %in% rownames(combinedLoadings), ] %>% t()
		testMatrix <- log2(testMatrix + 1) %>% as.matrix()
		
		predictedValues <- predict(model, ncomp = 2, newdata = testMatrix)
		
		tissueList1[[devTissue]] <- cor.test(predictedValues, posPSI)
		tissueList2[[devTissue]] <- data.frame(devTissue = devTissue, actualValue = posPSI, predictedValue = predictedValues)
	}
	cancerPredictionList[[cancerType]] <- tissueList1
	cancerPredictedList[[cancerType]] <- tissueList2
}

save(cancerPredictionList, cancerPredictedList, file = glue("predictionAndAccuracyAcorssTumors.Rda"))

predictionResults <- lapply(cancerPredictionList, function(x) lapply(x, function(y) c(correlation = y$estimate, pvalue = y$p.value))) %>% 
							lapply(., bind_rows, .id = "tissue") %>% bind_rows(., .id = "cancer")

dataForplot <- predictionResults %>% .[.$cancer == testCancers & .$tissue == names(testCancers), ]
dataForplot$tissue <- factor(dataForplot$tissue)
levels(dataForplot$tissue) <- gsub("Hindb", "B", levels(dataForplot$tissue))


barplot <- ggplot(data = dataForplot, aes(x = tissue, y = correlation.cor)) + theme_bw() + ylab("Correlation") + xlab("Cancer") +
				geom_bar(stat="identity", fill="white", col = "black") + 
				geom_text(aes(label = formatC(pvalue, format = "e", digits = 2),  size = 10, y = correlation.cor + 0.1), color = "black", position = position_dodge(0.9), size = 3.5)
								

ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 10),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 8),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
					
					
					
					
predictedResults <-  cancerPredictedList %>% lapply(., bind_rows, .id = "tissue") %>% bind_rows(., .id = "cancer")
dataForplot <- lapply(names(testCancers), function(x) predictedResults[predictedResults$cancer == testCancers[x] & predictedResults$tissue == x, ]) %>% do.call("rbind", .)
dataForplot$tissue <- factor(dataForplot$tissue)
levels(dataForplot$tissue) <- gsub("Hindb", "B", levels(dataForplot$tissue))
names(dataForplot)[5] <- "predictedValue"

ggSP <- ggplot(data = dataForplot, aes(x = actualValue, y = predictedValue)) + 
			theme_bw() + facet_wrap(~tissue, scale = "free", ncol = 1, strip.position="right") +
			#theme_bw() + facet_wrap(~tissue) +
			ylab("predicted inclusion") + xlab("actual inclusion") +
			geom_point(col = "mediumpurple2") + geom_smooth(method = "lm") + stat_cor(size = 2.5, label.sep='\n') + ggTheme
							

ggsave(file = glue("plotForPredictedSplicing.svg"), width =  2, height = 4) 