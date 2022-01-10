library(dplyr); library(magrittr); library(ensembldb); 
library(EnsDb.Hsapiens.v86); library(clusterProfiler)
library(gtools); library(reshape2); library(ggplot2)
library(clusterProfiler); library(gplots); library(glue)
orgDb <- "org.Hs.eg.db"; EnsdB <- EnsDb.Hsapiens.v86
library(TCGAbiolinks)

library(circlize)
library(ComplexHeatmap)


source('helperFunctions.R')
Tissue <- "Brain"; devTissue <- "Hindbrain"; cancerType <- "LGG"

Transcripts <- transcripts(EnsdB,
			          columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "tx_seq_start", "tx_seq_end"), #listColumns(EnsdB , "tx")),
			          #filter = TxBiotypeFilter("nonsense_mediated_decay"),
			          return.type = "DataFrame")
GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name))
#############load transcript expression files and collpase the values by gene expression level#############



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


keggPathwayFile <- download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")
pathwayToGene <- keggPathwayFile[[1]] %>% set_colnames(c("PathwayID", "EntrezID"))
pathwayToName <- keggPathwayFile[[2]] %>% set_colnames(c("PathwayID", "PathwayName"))


convertedIds <- bitr(pathwayToGene$EntrezID, "ENTREZID", "ENSEMBL", orgDb, drop = F)
matches <- match(pathwayToGene$EntrezID, convertedIds$ENTREZID)
pathwayToGene$GeneIDs <- convertedIds[matches, 'ENSEMBL']

matches <- match(pathwayToGene$GeneIDs, GenesDf$GeneID)
pathwayToGene$GeneName <- GenesDf[matches, 'GeneName']




load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))


desiredExons <- AssociationClusters[ ,'Exon'] %>% unique() %>% as.character()
AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
AssociationClusters <- AssociationClusters[AssociationClusters$ExonType != "Not defined", ]
AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	
embPos <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Positive", 'Exon'] %>% as.character()
embNeg <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Negative", 'Exon'] %>% as.character()



#colnames(psiPathwayCorrelation) <- rownames(devPSIfile)

fdrLevel <- 0.05
psiPathwayAssociation <- psiPathwayCorrelation[1:ncol(medianPathwayExpression), ] %>% t() %>% data.frame %>%
								set_rownames(rownames(devPSIfile)) %>% set_colnames(colnames(medianPathwayExpression))
								
psiPathwaySignificance <- psiPathwayCorrelation[(ncol(medianPathwayExpression) + 1):nrow(psiPathwayCorrelation), ]  %>% t() %>% data.frame %>%
								set_rownames(rownames(devPSIfile)) %>% set_colnames(colnames(medianPathwayExpression))
								
								
fdr <- apply(psiPathwaySignificance, 1, p.adjust, method = "fdr") %>% t() %>% data.frame


significanceMatrix <- apply(fdr, 2, function(x) ifelse(x <= fdrLevel, 1, 0))
psiPathwayAssociation <- psiPathwayAssociation * significanceMatrix
psiPathwayAssociation$Exon <- rownames(psiPathwayAssociation)
psiPathwayAssociation <- melt(psiPathwayAssociation)
names(psiPathwayAssociation) <- c("Exon", "Pathway", "Association")
psiPathwayAssociation <- psiPathwayAssociation[psiPathwayAssociation$Association != 0, ]

psiPathwayAssociation <- merge(psiPathwayAssociation, pathwayClusters, by.x = "Pathway", by.y = "pathwayId")
psiPathwayAssociation <- merge(psiPathwayAssociation, pathwayToName, by.x = "Pathway", by.y = "PathwayID")




correlatedPathwayDf <- psiPathwayAssociation[psiPathwayAssociation$clusterType == "Embryonic", ] %>% mutate(Correlation = ifelse(.$Association > 0, "Positive", "Negative"))
correlatedPathwayDf <- correlatedPathwayDf %>% group_by(., Exon, Correlation, clusterType) %>% summarize(Frequency = length(Correlation))
correlatedPathwayDf <- correlatedPathwayDf[ ,c(1,2,4)] %>% tidyr::spread(., key = Correlation, value = Frequency) %>% replace(., is.na(.), 0) %>% data.frame %>% set_rownames(.$Exon) %>% .[ ,-1] %>% data.matrix()


correlatedPathwayDf <- correlatedPathwayDf %>% .[rownames(.) %in% desiredExons, ]

exonGenesDf <- exonGeneNames(exonList = rownames(correlatedPathwayDf), GenesDf = GenesDf) %>% data.frame()


correlation <- list()
pvalue <- list()
for (index in 1:nrow(exonGenesDf)) {
	geneName <- exonGenesDf[index , 'GeneName'] %>% as.character()
	exonName <- exonGenesDf[index , 'Exon'] %>% as.character()
	
	psi <- devPSIfile[exonName, ] %>% unlist()
	expression <- devGeneTPM[geneName, ] %>% unlist()
	
	
	if (sum(is.na(expression)) > length(expression)/3) {
		correlation[[exonName]] <- NA
		pvalue[[exonName]] <- NA
	}else {
		res <- cor.test(psi, expression)
		correlation[[exonName]] <- res$estimate
		pvalue[[exonName]] <- res$p.value
	}
}

correlationDf <- data.frame(Exon = names(correlation), correlation = unlist(correlation), pvalue = unlist(pvalue)) %>% set_rownames(NULL)


#library(circlize)
#library(ComplexHeatmap)


my_palette <- colorRampPalette(c("violet", "white", "blue"))(n = 40)
col_fun1 = colorRamp2(c(0, 50, 100), c("violet", "white", "blue"))
col_fun2 = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))


lgd1 = Legend(title = "# correlated\npathways", col_fun = col_fun1)
lgd2 = Legend(title = "exon-gene\ncorrelation", col_fun = col_fun2)
lgd_list_vertical = packLegend(lgd1, lgd2, direction = "horizontal", row_gap = unit(2, "cm"))

set.seed(123)

indexesToplot <- sample(1:nrow(correlatedPathwayDf), 1000, replace = F)
subSetInclusion <- correlatedPathwayDf[indexesToplot, ]
subSetCorelation <- correlationDf[indexesToplot, ] %>% set_rownames(.$Exon) %>% subset(., select = correlation)


exonTypesDf <- rownames(subSetInclusion) %>% data.frame(exon = ., exonType = "Other Exons") %>% mutate(exonType = as.character(exonType))
exonTypesDf$exonType[exonTypesDf$exon %in% embPos] <- "Embryonic Positive"
exonTypesDf$exonType[exonTypesDf$exon %in% embNeg] <- "Embryonic Negative"
splitFactor <- factor(exonTypesDf$exonType)



#grid.draw(lgd1); grid.draw(lgd2)

#dev.copy(png, )
#dev.off()

tiff(filename = glue("exonCircularPlotFor{devTissue}.tiff"), width = 6, height = 6, units = "in", pointsize = 14,  bg = "white", res = 300, type = "quartz" , antialias = "none")
#bitmap(file = glue("exonCircularPlotFor{devTissue}.tiff"), type = "tiffg4",  width = 5, height = 5, units = "in", pointsize = 14,  bg = "white", res = 400)
circos.clear()
circos.heatmap(subSetInclusion, col = col_fun1, dend.side = "outside", split = splitFactor, show.sector.labels = TRUE)
circos.heatmap(subSetCorelation, col = col_fun2)
grid.draw(lgd_list_vertical)

dev.off()

save(correlationDf, correlatedPathwayDf, indexesToplot, file = glue("processedDataForCircularHeatmapIn{devTissue}.Rda"))

#heatmap.2(weightedCosine, dendrogram = 'row', Rowv = T, Colv = F, trace = 'none', scale = 'none', col = my_palette, density.info = 'none', keysize = 1, key.title = NA, key.xlab = "Cosine Similarity", xlab = "Developmental timeline", cexRow = 2, labRow = NA, labCol = NA, margins = c(2, 1), key.xtickfun = keyFunction, key.par=list(mgp=c(1, 0.5, 0)), lmat = lmat, lhei = lhei, lwid = lwid)



####make heatmap####



