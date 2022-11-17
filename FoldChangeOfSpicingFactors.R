#library(clusterProfiler)
library(dplyr); library(magrittr); library(ensembldb); library(EnsDb.Hsapiens.v86)
library(gtools); library(preprocessCore); library(reshape2); library(ggplot2)
library(glue); library(pls); library(TCGAbiolinks); library(data.table)
library(ggpubr); EnsdB <- EnsDb.Hsapiens.v86

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

		

usedTissues <- c("Hindbrain", "Liver", "Kidney")
tissueParameters <- list()
tissueParameters[['Hindbrain']] <- list(devTissue = "HindBrain", normalTissue = c("BrainCerebellum", "BrainCortex"), cancerType =  c("LGG", "GBM"), cancerStages = c("Overall"))
tissueParameters[['Liver']] <- list(devTissue = "Liver", normalTissue = "Liver", cancerType =  c("LIHC"), cancerStages = c("Overall", "Early", "Late"))
tissueParameters[['Kidney']] <- list(devTissue = "Kidney", normalTissue = "KidneyCortex", cancerType =  c("KIRP", "KIRC"), cancerStages = c("Overall", "Early", "Late"))


pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/"

foldChangeList <- list()
for (devTissue in usedTissues) {
	cancerType <- tissueParameters[[devTissue]]['cancerType'] %>% unlist %>% unname
	normalTissue <- tissueParameters[[devTissue]]['normalTissue'] %>% unlist %>% unname
	load(glue("plsrModelForMedianPosEmbSplicingIn{devTissue}.rda"))
	testGeneList <- combinedLoadings[combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0.0001, ] %>% rownames()
	#backgroundList <- rownames(combinedLoadings)
	#backgroundList <- backgroundList[!(backgroundList %in% testGeneList)]
	backgroundList <- combinedLoadings[combinedLoadings$coefficients <= 0.00, ] %>% rownames()
	
	
	
	devFile <- glue("{pathToExpressionFile}Normal/Kallisto/development/Expression/KallistoTPMvaluesForTranscriptExpressionIn{devTissue}")
	cancerExpression <- glue("{pathToExpressionFile}Cancer/Kallisto/TCGA/Expression/KallistoTPMvaluesForTranscriptExpressionIn{cancerType[1]}")
	if (devTissue == "Kidney") {
		cancerExpression <- glue("{pathToExpressionFile}Cancer/Kallisto/TCGA/Expression/KallistoTPMvaluesForTranscriptExpressionIn{cancerType[2]}")
	}
	normalExpression <- glue("{pathToExpressionFile}Cancer/Kallisto/Gtex/Expression/KallistoTPMvaluesForTranscriptExpressionIn{normalTissue[1]}")


	sampleNames <-  fread(paste(devFile), sep = "\t", header = F, check.names = F, skip = 0, nrow = 1) %>% as.character()
	devExpression <- read.table(devFile, sep = "\t", header = T, check.names = F, row.names = 1, colClasses=c("character", rep("numeric",length(sampleNames))))
	
	sampleNames <-  fread(paste(cancerExpression), sep = "\t", header = F, check.names = F, skip = 0, nrow = 1) %>% as.character()
	cancerExpression <- read.table(cancerExpression, sep = "\t", header = T, check.names = F, row.names = 1, colClasses=c("character", rep("numeric",length(sampleNames))))
		
		
	sampleNames <-  fread(paste(normalExpression), sep = "\t", header = F, check.names = F, skip = 0, nrow = 1) %>% as.character()
	normalExpression <- read.table(normalExpression, sep = "\t", header = T, check.names = F, row.names = 1, colClasses=c("character", rep("numeric",length(sampleNames))))

	names(cancerExpression) <- gsub("\\.", "-", names(cancerExpression))
	names(normalExpression) <- gsub("\\.", "-", names(normalExpression))

	cancerSamples <- TCGAquery_SampleTypes(names(cancerExpression), typesample = c("TP"))
	cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples]

	commonTranscripts <- rownames(cancerExpression)[rownames(cancerExpression) %in% rownames(normalExpression)] %>% .[. %in% rownames(devExpression)]
						
	cancerExpression <- cancerExpression[rownames(cancerExpression) %in% commonTranscripts, ] %>% .[order(rownames(.)), ]
	normalExpression <- normalExpression[rownames(normalExpression) %in% commonTranscripts, ] %>% .[order(rownames(.)), ]
	devExpression <- devExpression[rownames(devExpression) %in% commonTranscripts, ] %>% .[order(rownames(.)), ]

	combinedFile <- cbind(cancerExpression, normalExpression)
	geneLevelTPM <- processEpxressionFile(dataFrame = combinedFile, Transcripts = Transcripts, replicates = F)

	cancerGeneTPM <- geneLevelTPM[ ,names(geneLevelTPM) %in% names(cancerExpression)]
	normalGeneTPM <- geneLevelTPM[ ,names(geneLevelTPM) %in% names(normalExpression)]

	
	normalMedian <- apply(normalGeneTPM, 1, median, na.rm = T)
	cancerMedian <- apply(cancerGeneTPM, 1, median, na.rm = T)

	#foldChange <- log2((cancerMedian / normalMedian) + 1)
	foldChange <- log2(cancerMedian / normalMedian)
 

	criticalFoldChange <- foldChange[names(foldChange) %in% testGeneList]
	OtherFoldchange <-  foldChange[names(foldChange) %in% backgroundList]

	#combinedLoadings$factorType <- ifelse(combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0, "Positive Regulators", "Other Factors")
	#combinedLoadings$factorType[combinedLoadings$coefficients < 0 & combinedLoadings$FDR < 0.05] <- "Negative Regulators"
	
	#combinedLoadings$factorType <- ifelse(combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0, "Critical", "Not Critical")
	
	combinedLoadings <- merge(combinedLoadings, foldChange, by = 0)
	names(combinedLoadings)[names(combinedLoadings) == "y"] <- "cancerFoldChange"
	
	

	devGeneTPM <- processEpxressionFile(dataFrame = devExpression, Transcripts = Transcripts, replicates = T)
	devGeneTPM <- devGeneTPM[ ,mixedsort(names(devGeneTPM))]

	stage <- gsub("\\w+_|days\\w?", "", names(devGeneTPM)) %>% as.numeric()
	preNatal <- stage < 300
	postNatal <- stage > 300
	
	preNatalMedian <- apply(devGeneTPM[ ,preNatal], 1, median, na.rm = T)
	postNatalMedian <- apply(devGeneTPM[ ,postNatal], 1, median, na.rm = T)
	foldChange <- log2(preNatalMedian / postNatalMedian)
	combinedLoadings <- merge(combinedLoadings, foldChange, by.x = "Row.names",  by.y = 0)
	names(combinedLoadings)[names(combinedLoadings) == "y"] <- "developmentFoldChange"
	
	combinedLoadings$factorType <- ifelse(combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0, "Critical", "Not Critical")
	foldChangeList[[devTissue]] <- combinedLoadings
}
save(foldChangeList, file = "foldChangeListOfSplicingFactors.Rda")


dataForplot <- bind_rows(foldChangeList, .id = "Cancer")
dataForplot <- reshape2::melt(dataForplot, id.vars=c("factorType", "Cancer"), measure.vars = c("cancerFoldChange", "developmentFoldChange"))
colnames(dataForplot) <- c("factorType", "Cancer", "FoldChangeType", "FoldChange")


ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 15, vjust = 1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 10),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 10),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	

dataForplot$FoldChangeType <- factor(dataForplot$FoldChangeType, levels = c("developmentFoldChange", "cancerFoldChange"))
levels(dataForplot$FoldChangeType) <- c("Pre- Vs. Post-Natal", "Cancer Vs Normal")
ggplot(data = dataForplot, aes(x = factorType, y = FoldChange)) + facet_grid(Cancer~FoldChangeType) + theme_bw() +  ggTheme +
			geom_boxplot(outlier.shape = NA) + stat_compare_means(label = "p.format", label.x = 1.4, label.y = 3, size = 4) + coord_cartesian(ylim = c(-4,4))

ggsave(file = glue("plotForFoldChangeOfSplicingFactors.svg"), width =  3, height = 4) 

