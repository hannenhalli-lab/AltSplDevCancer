#library(clusterProfiler)
library(dplyr); library(magrittr); library(ensembldb); library(EnsDb.Hsapiens.v86)
library(gtools); library(preprocessCore); library(reshape2); library(ggplot2)
library(clusterProfiler); library(msigdbr); library(gplots); library(FactoMineR); library(factoextra)
library(TCGAbiolinks); library(glue); orgDb <- "org.Hs.eg.db"; EnsdB <- EnsDb.Hsapiens.v86

Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),

GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))

pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/Normal/Kallisto/"
tissueCancerFile <- read.table("../Cancer/Kallisto/cancerAndTissuePairs", sep = "\t", header = T)




processPSIfile <- function(psiFile) {
	data.frame(t(psiFile)) %>% 
	mutate(Stage = gsub("_Rep\\d+", "", names(psiFile))) %>% 
	aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
	set_rownames(.$Stage) %>% 
	.[,-1] %>% t() %>% data.frame() %>%
	set_rownames(rownames(psiFile))
}

exonGeneNames <- function(exonList, GenesDf) {
	exonGenes <-  exonList %>% as.character() %>% gsub("\\..+", "", .)
	indexes <- match(exonGenes, GenesDf[ ,'GeneID'])
	exonGenes <- GenesDf[indexes, c('GeneName')] %>% as.character()
	geneType <- GenesDf[indexes, c('GeneType')] %>% as.character()
	return(cbind(exonGenes, geneType))
}

commonExonGenes <- function(list1, list2) {
	list1 <- as.character(list1) %>% unique(); list2 <- as.character(list2) %>% unique()
	commonGenes <- list1[list1 %in% list2]
	return(commonGenes)
}

tissueTypes <- tissueCancerFile[ ,'devTissue'] %>% as.character() %>% unique()

medianSplicingDf <- NULL
for (fileIndex in 1:length(tissueTypes)) {
	devTissue <- tissueTypes[fileIndex]
	
	load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))
	
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	#AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType != "Not defined", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	
	embPos <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Positive", 'Exon'] %>% as.character()
	embNeg <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Negative", 'Exon'] %>% as.character()
	
	posGeneNames <- exonGeneNames(exonList = embPos, GenesDf) %>% data.frame()
	negGeneNames <- exonGeneNames(exonList = embNeg, GenesDf) %>% data.frame()
	
	commonGenes <- commonExonGenes(list1 = posGeneNames$exonGenes, list2 = negGeneNames$exonGenes)
	
	commonIndexes <- posGeneNames$exonGenes %in% commonGenes
	embPos <- embPos[commonIndexes]
	
	commonIndexes <- negGeneNames$exonGenes %in% commonGenes
	embNeg <- embNeg[commonIndexes]
		
	devFile <- glue("{pathToExpressionFile}development/PSIvalues/SuppaPSIvaluesUsingKallistoIn{devTissue}.psi")
	devPSIfile <- read.table(devFile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)
	
	NAsums <- apply(devPSIfile,1, function(x) sum(is.na(x)))
	NAthreshold <- ceiling(ncol(devPSIfile) / 2)
	indexes <- which(NAsums < NAthreshold)
	devPSIfile <- devPSIfile[indexes, ]

	devPSIfile <- processPSIfile(psiFile = devPSIfile)
	devPSIfile <- devPSIfile[, mixedsort(names(devPSIfile))]			
	
	posPSI <- devPSIfile[rownames(devPSIfile) %in% embPos, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
	negPSI <- devPSIfile[rownames(devPSIfile) %in% embNeg, ] %>% apply(., 2, mean, na.rm = T, na.action = na.pass)
	
	#medianPSI <- apply(devPSIfile, 2, mean, na.rm = T)
	
	ggTemp <- data.frame(sample = names(posPSI), embryonicPos = posPSI, embryonicNeg = negPSI, tissue = devTissue)
	medianSplicingDf <- rbind(medianSplicingDf, ggTemp)
}
rownames(medianSplicingDf) <- NULL
medianSplicingDf <- melt(medianSplicingDf, id.vars = c("sample", "tissue")) %>% set_colnames(c("sample", "tissue", "exonType", "InclusionLevel"))


ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 45, vjust =1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	legend.title = element_blank(), legend.background = element_blank(),
	axis.title = element_text(size = 12, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 1))

lineplotFunction <- function(data) {
	ggplot(data = data, aes(x = sample, y = InclusionLevel, group = exonType)) + theme_bw() + ggTheme + theme(legend.position = c(0.5, 0.8)) + 
			geom_point(aes(color = exonType)) + geom_line(aes(color = exonType)) + facet_wrap(~tissue, ncol = 1, scale = "free") 
							
}

boxplotFunction <- function(data, tissue) {
	ggplot(data = data, aes(x = stage, y = InclusionLevel)) + theme_bw() + geom_boxplot() + facet_wrap(~exonType, scale = "free") + ggtitle(paste(tissue)) +
		ggTheme + theme(axis.text.x = element_text(size = 10, angle = 20, vjust = 1, hjust = 1))										
}					

plotTissues <- tissueTypes[1:3]	
for (Tissue in plotTissues) {
	subData <- medianSplicingDf[medianSplicingDf$tissue == Tissue, ]
	subData$sample <- gsub(glue("{Tissue}_"), "", subData$sample)
	factorOrder <- unique(subData$sample) %>% as.character() 
	subData$sample <- factor(subData$sample, levels = factorOrder)
	subData$stage <- gsub("days", "", subData$sample) %>% as.character() %>% as.numeric()
	subData$stage <- ifelse(subData$stage  < 200, "Pre Natal", "Post Natal")
	subData <- melt(subData) %>% set_names(c("sample", "tissue", "stage", "exonType", "InclusionLevel"))
	#plotData <- subData[subData$exonType == "embryonicPos", ]
	plotData <- subData
	plotData$stage <- factor(plotData$stage, levels = c("Pre Natal", "Post Natal"))
	boxplotFunction(data = plotData, tissue = Tissue)
	ggsave(glue('developmentalSplicingBoxplotFor{Tissue}.svg'), width =  3.2, height = 3.2)
}
#Tissue = "Kidney"


