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
	mutate(depMapID = gsub("-", "\\.", depMapID)) %>% subset(. , select = c(depMapID, Name, Doubling.Time.from.Vendor, Doubling.Time.Calculated.hrs))





devTissues <- c("Hindbrain" = "Brain Cancer", Kidney = "Kidney Cancer", Liver = "Liver Cancer")
#desiredTypes <- list(Hindbrain = c("Glioblastoma", "Glioma"), Kidney = "all", Liver = "all")
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

	subData <- merge(subSamples, medianPosSplicing, by.x = "DepMap_ID", by.y = 0) %>% merge(., medianNegSplicing, by.x = "DepMap_ID", by.y = 0) %>% 
				set_names(c("DepMap_ID", "Subtype", "stage", "primary_disease", "vendorTime", "calculatedTime", "embryonic.Positive", "embryonic.Negative"))
	
	
	
	dataForplotList[[cancer]] <- subData
}
save(dataForplotList, file = "ccleProliferationAndEmbryonicSplicing.Rda")

ggTheme <- theme(axis.text.x = element_text(size = 8, angle = 20, vjust = 1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 10, face = "bold"),
	#strip.background = element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	plot.margin = unit(c(0, 0, 0, 0), "cm"),
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
	
	
dataForplot <- bind_rows(dataForplotList, .id = "Cancer")
dataForplot$Cancer <- factor(dataForplot$Cancer)

dataForplot <- dataForplot %>% subset(., select = c(Cancer, DepMap_ID, calculatedTime, embryonic.Positive, embryonic.Negative)) %>%
	 			reshape2::melt(., id.vars = c("DepMap_ID", "calculatedTime", "Cancer")) %>% set_names(c("DepMap_ID", "calculatedTime", "Cancer", "exonType", "inclusionLevel"))

levels(dataForplot$exonType) <- c("EP", "EN")

dataForplot$Cancer <- gsub(" Cancer", "", dataForplot$Cancer)

p1 <- dataForplot %>% subset(., exonType == "EP") %>%  
 		ggplot(data = ., aes(x = inclusionLevel, y = calculatedTime)) + facet_wrap(~Cancer, scale = "free", ncol = 1) + 
		ylab("Doubling time") + xlab("Inclusion Level") +
		geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman", label.x.npc = 0, label.y.npc = 0.8, label.sep='\n') + 
		theme(strip.background = element_blank(), strip.text = element_blank()) + ggTheme


p2 <- dataForplot %>% subset(., Cancer != " Cancer" & exonType == "EN") %>%  
 		ggplot(data = ., aes(x = inclusionLevel, y = calculatedTime)) + 
		facet_wrap(~Cancer, scale = "free", strip.position = c("right"), ncol = 1) + 
		ylab("")  + xlab("Inclusion Level") +
		geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman", label.x.npc = 0, label.y.npc = 0.8, label.sep='\n') + 
		theme(strip.background = element_blank(), strip.text = element_blank()) + ggTheme
		#theme(strip.background = element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12)) + ggTheme



#library(cowplot)
ggBP <- plot_grid(p1, p2)
ggsave(ggBP, file = glue("ccleProliferationRatesAndInclusionLevel.svg"), width =  3.5, height = 5.5) 



