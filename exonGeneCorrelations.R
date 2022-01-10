library(TCGAbiolinks); library(glue); 
library(dplyr); library(data.table);
library(maftools); library(EnsDb.Hsapiens.v86)
library(GenomicFeatures); library(magrittr)
library(TCGAbiolinks); library(OneR)
EnsdB <- EnsDb.Hsapiens.v86
library(ggplot2); library(ggpubr); library(gtools)
#library(reshape2)


getEvents <- function (AssociationClusters, ExonType, Association) {
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	ExonType <- paste(ExonType, Association, sep = ".")
	events <- AssociationClusters[AssociationClusters$ExonType == ExonType, 'Exon'] %>% as.character()
	return(events)
}


source('helperFunctions.R')

Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),

GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))



devTissue <- "Kidney"
load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))

embPos <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Positive")
embNeg <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Negative")

exonDf <- rbind(data.frame(Exon = embPos, exonType = "EmbryonicPositive"), data.frame(Exon = embNeg, exonType = "EmbryonicNegative"))

pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/"

devFile <- glue("{pathToExpressionFile}Normal/Kallisto/development/PSIvalues/SuppaPSIvaluesUsingKallistoIn{devTissue}.psi")
devPSIfile <- read.table(devFile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)

devPSIfile <- naFilter(dataFrame = devPSIfile, cutoff = 0.5)
devPSIfile <- processPSIfile(psiFile = devPSIfile)
devPSIfile <- devPSIfile[, mixedsort(names(devPSIfile))]



devExpressionFile <- glue("{pathToExpressionFile}Normal/Kallisto/development/Expression/KallistoTPMvaluesForTranscriptExpressionIn{devTissue}")
devExpression <- read.table(paste(devExpressionFile), sep = "\t", header = T, check.names = F)

devGeneTPM <- processEpxressionFile(dataFrame = devExpression, Transcripts = Transcripts, replicates = T)
devGeneTPM <- devGeneTPM[ ,mixedsort(names(devGeneTPM))]


splicingSubData <- devPSIfile %>% .[rownames(.) %in% exonDf[['Exon']], ]
exonGenesDf <- exonGeneNames(exonList = rownames(splicingSubData), GenesDf = GenesDf)
exonGenesDf <- merge(exonGenesDf, exonDf, by = "Exon")


correlation <- list()
pvalue <- list()
for (index in 1:nrow(exonGenesDf)) {
	geneName <- exonGenesDf[index , 'GeneName'] %>% as.character()
	exonName <- exonGenesDf[index , 'Exon'] %>% as.character()
	
	psi <- splicingSubData[exonName, ] %>% unlist()
	expression <- devGeneTPM[geneName, ] %>% unlist()
	
	res <- cor.test(psi, expression)
	correlation[[exonName]] <- res$estimate
	pvalue[[exonName]] <- res$p.value
}


correlationDf <- data.frame(Exon = names(correlation), correlation = unlist(correlation), pvalue = unlist(pvalue)) %>% set_rownames(NULL)
correlationDf$correlationType <- ifelse(correlationDf$correlation > 0, "positive", "negative")
correlationDf$correlationType[correlationDf$pvalue > 0.1] = "uncorrelated"
correlationDf <- merge(correlationDf, exonGenesDf, by = "Exon")

save(correlationDf, file = glue("exonGeneCorrelationForEmbryonicExonsIn{devTissue}.Rda"))


#correlationDf$correlationType <- factor(correlationDf$correlationType, levels = c("negative", "uncorrelated", "positive"))
#boxplot(correlation ~ correlationType, data = correlationDf)