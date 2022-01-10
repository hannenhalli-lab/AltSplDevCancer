###############use knockTF data#################
#library(clusterProfiler)
library(dplyr); library(magrittr); library(gplots); library(glue)
library(gtools); library(reshape2); library(ggplot2);
library(metap) ## for fisher combination

#load('plsrModelForMedianPosEmbSplicingInHindbrain.rda')

hypergeometricFunction <- function(list1, list2, universe) {
	
	list1 <- list1[list1 %in% universe]
	list2 <- list2[list2 %in% universe]
	universe <- length(universe)
	
	#c1 <- subData2[subData2$ExonType == exonType, 'Count']
	c1 <- sum(list1 %in% list2)
	#c2 <- backgroundInCancer[eventType] - c1
	c2 <- length(list1) - c1
	
	#c3 <- backGroundInPathways[exonType] - c1
	c3 <- length(list2) - c1
	
	#c4 <- universe - (c3 + backgroundInCancer[eventType])
	c4 <- universe - (c3 + length(list1))
	
	#hyperP <- phyper(c1 - 1, backgroundInCancer[eventType], universe - backgroundInCancer[eventType], backGroundInPathways[exonType], lower.tail= FALSE)
	hyperP <- phyper(c1 - 1, length(list1), universe - length(list1), length(list2), lower.tail= FALSE)
	
	return(c(sizeOfTargets = length(list1), sizeOfInput = length(list2), overlap = c1, pvalue = hyperP))
}

fisherTestFunction <- function(targets, testGeneList, backgroundList) {
	testOverlap <- table(targets %in% testGeneList)
	testOverlap <- testOverlap['TRUE']
	if (is.na(testOverlap)) {
		testOverlap <- 0
	}
	testNonOverlap <- length(testGeneList) - testOverlap
	

	backOverlap <- table(targets %in% backgroundList)
	backOverlap <- backOverlap['TRUE']
	if (is.na(backOverlap)) {
		backOverlap <- 0
	}
	backNonOverlap <- length(backgroundList) - backOverlap
	
	fisherMat <- matrix(c(testOverlap, testNonOverlap, backOverlap, backNonOverlap), nrow = 2)
	testRes <- fisher.test(fisherMat)
	oddsRatio <- testRes$estimate
	pvalue <- testRes$p.value

	testOverlap <- testOverlap / length(testGeneList)
	backOverlap <- backOverlap / length(backgroundList)

	tempDf <- cbind(testOverlap = testOverlap, backgroundOverlap = backOverlap, oddsRatio = oddsRatio, pvalue = pvalue) %>% data.frame()
	return(tempDf)
}


devTissues <- c("Hindbrain", "Kidney", "Liver")

#colClasses <- c(rep("character", 3), rep("numeric", 7))
colClasses <- c(rep("character", 10))
knockTF <- read.table("knockTFdata/differential expression of genes in all datasets.txt", sep = "\t", header = T)#, colClasses = colClasses)
#knockTF <- knockTF %>% varhandle::unfactor(knockTF)
datasetInfo <- read.csv("knockTFdata/KnockTF-Browse.csv")
uniqueTargets <- knockTF$Gene %>% unique() %>% as.character()
allTFs <-  knockTF$TF %>% unique() %>% as.character()
allTargets <- knockTF$Gene %>% unique() %>% as.character()


load('topTranscriptionFactorList.Rda')


tfEnrichmentList <- list()
for (devTissue in devTissues) {
	load(glue('plsrModelForMedianPosEmbSplicingIn{devTissue}.rda'))
	tissue <- gsub("Hindb", "B", devTissue)

	testGeneList <- combinedLoadings[combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0.0001, ] %>% rownames() #%>% GeneID2entrez(.)
	#backgroundList <- rownames(combinedLoadings) #%>% GeneID2entrez(.)
	#backgroundList <- backgroundList[!(backgroundList %in% testGeneList)]
	backgroundList <- combinedLoadings[combinedLoadings$coefficients < 0.0, ] %>% rownames()
	

	####select again top genes####
	if (devTissue == "Kidney") {
		testGeneList <- combinedLoadings %>% .[.$FDR <= 0.05, ] %>% .[order(.$coefficients, decreasing = T), ] %>% rownames %>% .[1:45] 
	} else {
		#testGeneList <- combinedLoadings %>% .[.$FDR <= 0.05, ] %>% .[order(.$coefficients, decreasing = T), ] %>% rownames %>% .[1:100]
	}

	#testGeneList <- testGeneList[testGeneList %in% allTargets]
	#backgroundList <- backgroundList[backgroundList %in% allTargets]
	
	topTFlist <- topTFs[[tissue]]
	enrichmentList <- list()
	for (index in 1:length(topTFlist)) {
		TF <- topTFlist[index]
		subData <- knockTF[knockTF$TF == TF, ]
		length(unique(subData$Sample_ID))
		subData <- subData[subData$up_down == 2, ]
		subData$Sample_ID <- factor(subData$Sample_ID, levels = unique(subData$Sample_ID))
	
		if (nrow(subData) > 5) {
			#gg1 <- subData %>% split(., .$Sample_ID) %>% lapply(., function(x) length(unique(x$Gene)))
			fisherResult <- subData %>% split(., .$Sample_ID) %>% lapply(., function(x) 
					fisherTestFunction(targets = x$Gene, testGeneList = testGeneList, backgroundList = backgroundList))
				
			fisherResultDf <- bind_rows(fisherResult, .id = "SampleID")
			enrichmentList[[TF]] <- fisherResultDf
		}
	}
	enrichmentDf <- bind_rows(enrichmentList, .id = "TF")
	enrichmentDf$FDR <- p.adjust(enrichmentDf$pvalue, method = "fdr")
	tfEnrichmentList[[tissue]] <- enrichmentDf
}

tfEnrichmentDf <-  bind_rows(tfEnrichmentList, .id = "Tissue")

save(tfEnrichmentDf, file = glue('knockTFenrichmentForTFs-finalRun.Rda'))

tfEnrichmentDf$significance <- ifelse(tfEnrichmentDf$FDR < 0.3 , "sig", "nsig")
tfEnrichmentDf$label <- "NoLabel"
tfEnrichmentDf$label[tfEnrichmentDf$FDR < 0.3 & tfEnrichmentDf$oddsRatio > 1] <- "Enriched"


fileName <- "knockTFresultsForTissues.svg"
ggTheme <- theme(axis.text.x = element_text(size = 14), #face = "bold"),
                 axis.text.y = element_text(size = 14, face = "bold"), plot.title = element_text(size = 14, hjust = 0.5),
                 legend.title = element_blank(), #legend.text = element_text(size = 30), legend.key.size = unit(2, 'cm'), 
				 legend.position = 'none',
                 axis.title = element_text(size = 14, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 1))

plot1 <- ggplot(data = tfEnrichmentDf, aes(x = Tissue, y = log2(oddsRatio)), color = significance) + geom_point() + geom_vline(xintercept = c(1, 2, 3)) + theme_bw() + ggTheme
plot1 <- plot1 + geom_label_repel(data = subset(tfEnrichmentDf, label == "Enriched"), aes(label = TF, color = label), 
    size = 2.8, fontface = 'bold', box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), nudge_y = 0.5, max.overlaps = 30) +
    coord_cartesian(ylim = c(-5,5)) + ggtitle("knockTF")

#plot1 <- plot1 + guides(color = guide_legend(override.aes = list(size = 12)))
ggsave(plot1, file = fileName, height = 6, width = 6)

		


