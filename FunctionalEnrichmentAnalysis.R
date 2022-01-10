library(clusterProfiler)
library(org.Hs.eg.db)
orgDb <- "org.Hs.eg.db"
library(dplyr); library(magrittr);
library(glue); library(ggplot2)

#####load PCA results######
#load('pcaResultsForLM22deconvolutionFromSaleem.Rda')


idConversion <- function(dataForTest, from, to, orgDb, removeNA) {
	geneList <- dataForTest #names(dataForTest)
	idsDf <- bitr(geneList, fromType = from, toType = to, OrgDb = orgDb, drop = F)
	idsDf <- idsDf[!duplicated(idsDf[,c('SYMBOL')]),] ###randomly retain unique ids from one to multiple comversion cases
	geneList <- idsDf[[to]]
	if (removeNA) {
		geneList <- geneList[!(is.na(geneList))]
	}
	return(geneList)
}


getPathwayGenes <- function(pathwayClusters, pathwayToGene, clusterType, ranked, memberShipCutoff) {
	names(pathwayClusters) <- c("Pathway", "Cluster", "ClusterType")
	names(pathwayToGene) <- c("Pathway", "EntrezID")
	geneClusters <- merge(pathwayClusters, pathwayToGene, by = "Pathway", all.x = T)
	
	if (!ranked) {
		geneList <- geneClusters[geneClusters$ClusterType == clusterType, ] %>%  
						group_by(., EntrezID) %>% summarise(memberShipSize = length(Pathway)) #%>% 
						#subset(., memberShipSize >= memberShipCutoff, select = EntrezID) %>% unique %>% 
						#data.frame %>% unlist	
	}		
	
	if (ranked) {
		geneList <- geneClusters[geneClusters$ClusterType == clusterType, ] %>%  
						group_by(., EntrezID) %>% summarise(memberShipSize = length(Pathway)) %>% 
						tibble::deframe() %>% .[order(., decreasing = T)]
	}
	return(geneList)		
}

keggPathwayFile <- download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")
pathwayToGene <- keggPathwayFile[[1]] %>% set_colnames(c("PathwayID", "EntrezID"))
pathwayToName <- keggPathwayFile[[2]] %>% set_colnames(c("PathwayID", "PathwayName"))
universe <- pathwayToGene$EntrezID %>% unique


devTissue <- "Kidney"
load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))




embryonic <- getPathwayGenes(pathwayClusters = pathwayClusters, pathwayToGene = pathwayToGene, clusterType = "Embryonic", ranked = F, memberShipCutoff = 5)
nonEmbryonic <-  getPathwayGenes(pathwayClusters = pathwayClusters, pathwayToGene = pathwayToGene, clusterType = "Adult", ranked = F, memberShipCutoff = 5)
#nonEmbryonic <- geneClusters[geneClusters$clusterType == "Adult", 'EntrezID'] %>% unique

combinedDf <- merge(embryonic, nonEmbryonic, by = "EntrezID", all = T) %>% set_colnames(c("EntrezID", "embryonicSize", "adultSize"))
combinedDf$foldChange <- combinedDf$embryonicSize / combinedDf$adultSize 

embryonicIndexes <- which(combinedDf$embryonicSize > 5 & combinedDf$foldChange >= 2)
nonEmbryonicIndexes <- which(combinedDf$adultSize > 5 & combinedDf$foldChange <= 0.5)


embryonicGenes <- combinedDf[embryonicIndexes, 'EntrezID']
nonEmbryonicGenes <- combinedDf[nonEmbryonicIndexes, 'EntrezID']


#commonGenes <- embryonic[embryonic %in% nonEmbryonic]
#embryonicGenes <- embryonic[!(embryonic %in% commonGenes)]
#nonEmbryonicGenes <- nonEmbryonic[!(nonEmbryonic %in% commonGenes)]


goEnrihmentResult <- list()

#interestingKeywords <- c("epithelial", "cytoskeleton", "proliferation", " stem ", "progenitor", "cell cycle", "mitosis", "angiogenesis", "adhesion", "development", "differentiation", "tissue", "signal", "cancer", "tumor", "apoptosis", "notch", "wnt")
#interestingKeywords <- paste(interestingKeywords,collapse="|")

enrichment <- enrichGO(embryonicGenes,  ont ="BP",  keyType = "ENTREZID", OrgDb = orgDb, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, universe = universe)
level4 <- gofilter(enrichment, level = 4)
sementic <- simplify(level4, cutoff = 0.7)

#dropTerms <- sementic %>% data.frame %>% .[grep(interestingKeywords, .$Description, ignore.case = T, invert = T), 'ID']
#forPlot <- dropGO(sementic, term = dropTerms)
goEnrihmentResult[['embryonic']] <- c(enrichment = enrichment, level4 = level4, sementic = sementic)


#dotplot(forPlot, color="qvalue", showCategory = 50) + ggtitle("embryonic") + theme(axis.text.y = element_text(size = 8, face = "bold"))


enrichment <- enrichGO(nonEmbryonicGenes,  ont ="BP",  keyType = "ENTREZID", OrgDb = orgDb, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, universe = universe)
level4 <- gofilter(enrichment, level = 4)
sementic <- simplify(level4, cutoff = 0.7)

#dropTerms <- sementic %>% data.frame %>% .[grep(interestingKeywords, .$Description, ignore.case = T, invert = T), 'ID']
#forPlot <- dropGO(sementic, term = dropTerms)

goEnrihmentResult[['nonEmbryonic']] <- c(enrichment = enrichment, level4 = level4, sementic = sementic)







#interestingKeywords <- c("epithelial", "cytoskeleton", "proliferation", " stem ", "progenitor", "cell cycle", "mitosis", "angiogenesis", "adhesion", "development", "differentiation", "tissue", "signal", "cancer", "tumor", "apoptosis", "notch", "wnt")


keggPathwayResults <- list()

enrichment1 <- enrichKEGG(embryonicGenes, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, universe = universe)
#dropTerms <- enrichment1 %>% data.frame %>% .[grep(interestingKeywords, .$Description, ignore.case = T, invert = T), 'ID']
#forPlot <- dropGO(enrichment1, term = dropTerms)

enrichment2 <- enrichKEGG(nonEmbryonicGenes, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, universe = universe)
#dropTerms <- enrichment2 %>% data.frame %>% .[grep(interestingKeywords, .$Description, ignore.case = T, invert = T), 'ID']
#forPlot <- dropGO(enrichment2, term = dropTerms)



#commonTerms <- enrichment1[enrichment1$Description %in% enrichment2$Description, 'ID']
#forPlot1 <- dropGO(enrichment1, term = commonTerms)
#forPlot2 <- dropGO(enrichment2, term = commonTerms)



keggPathwayResults[['embryonic']] <- c(enrichment = enrichment1)
keggPathwayResults[['nonEmbryonic']] <- c(enrichment = enrichment2)
save(goEnrihmentResult, keggPathwayResults, file = glue('functionalTermsOfembryonicPathwaysIn{devTissue}.Rda'))



dotplot(forPlot1, color="qvalue", showCategory = 80) +  theme(axis.text.y = element_text(size = 8, face = "bold"))
embryonic <- getPathwayGenes(pathwayClusters = pathwayClusters, pathwayToGene = pathwayToGene, clusterType = "Embryonic", ranked = T, memberShipCutoff = 10)
nonEmbryonic <-  getPathwayGenes(pathwayClusters = pathwayClusters, pathwayToGene = pathwayToGene, clusterType = "Adult", ranked = T, memberShipCutoff = 10)
ego3 <- gseGO(geneList = embryonic[embryonic > 1],  OrgDb = orgDb, ont = "BP", nPerm = 500, minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH")
bp1 <- gofilter(ego3, level = 3)
bp1 <- simplify(bp1, cutoff=0.7)











################################################
devTissue <- "Hindbrain"
load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))


AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
AssociationClusters <- AssociationClusters[AssociationClusters$ExonType != "Not defined", ]
AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	
embPos <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Positive", 'geneName'] %>% as.character() %>% unique()
embNeg <- AssociationClusters[AssociationClusters$ExonType == "Embryonic.Negative", 'geneName'] %>% as.character() %>% unique()
commonGenes <- embPos[embPos %in% embNeg]

#embPos <- embPos[!(embPos %in% commonGenes)]
#embNeg <- embNeg[!(embNeg %in% commonGenes)]



embPos <- idConversion(dataForTest = embPos, from = "SYMBOL", to = "ENTREZID", orgDb = org.Hs.eg.db, removeNA = T) %>% unique()
embNeg <- idConversion(dataForTest = embNeg, from = "SYMBOL", to = "ENTREZID", orgDb = org.Hs.eg.db, removeNA = T) %>% unique()



goEnrihmentResult <- list() 
keggPathwayResults <- list()
ggGenes1 <- c(embPos, embNeg) %>% unique

enrichment <- enrichGO(ggGenes1,  ont ="BP",  keyType = "ENTREZID", OrgDb = orgDb, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
level4 <- gofilter(enrichment, level = 4)
sementic <- simplify(level4, cutoff=0.7)
goEnrihmentResult[['devExons']] <- c(enrichment = enrichment, level4 = level4, sementic = sementic)
#dotplot(bp1, showCategory = 70)

#dotplot(bp1, color="qvalue", showCategory = 50) + ggtitle(glue(test)) + theme(axis.text.y = element_text(size = 8, face = "bold"))

####################################################################################################################################
enrichment <- enrichKEGG(ggGenes1, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
keggPathwayResults[['devExons']] <- list(enrichment = enrichment)


save(goEnrihmentResult, keggPathwayResults, file = glue('functionalTermsOfembryonicExons{devTissue}.Rda'))


#dotplot(enrichment, color="qvalue", showCategory = 50) + ggtitle(glue(test)) + theme(axis.text.y = element_text(size = 8, face = "bold"))



enrichment <- enrichGO(ggGenes1,  ont ="MF",  keyType = "ENTREZID", OrgDb = orgDb, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
level4 <- gofilter(enrichment, level = 4)
sementic <- simplify(level4, cutoff=0.7)
goEnrihmentResult[['devExons']] <- c(enrichment = enrichment, level4 = level4, sementic = sementic)

save(goEnrihmentResult, keggPathwayResults, file = glue('MolecularFunctionsOfembryonicExons{devTissue}.Rda'))



