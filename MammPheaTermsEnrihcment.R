library(dplyr); library(magrittr); library(glue)
library(reshape2); library(ggplot2); library(data.table)
library(tidyverse); library(ontologyIndex)


#termDefinitions <- read.table("VOC_MammalianPhenotype.rpt.txt", sep = "\t", header = F)
#termMappings <- fread("HMD_HumanPhenotype.rpt.txt") %>% set_colnames(c("geneName", "enterzID", "mouseGene", "MGI", "terms", "blank")) %>% .[ ,1:5]
#termMappings <- termMappings %>% mutate(terms = strsplit(terms, ",")) %>% unnest(cols = terms) %>% mutate(terms = gsub(" +", "", .$terms))

termMappings <- read.table("gene_attribute_edges_simplified.txt", sep = "\t", header = T, skip = 1)
allGenes <- termMappings$GeneSym %>% unique %>% as.character()

ontologyFile <- get_ontology("MPheno_OBO.ontology")
termDefinitions <- ontologyFile %>% as.data.frame() %>% select(c(id, name, obsolete, children))


#embryonicPhenotypes <- c("embryo", "prewean", "lethal") %>% paste(., collapse = "|")
#embryonicTerms <- termDefinitions[grep(embryonicPhenotypes, termDefinitions$name), ]

listOfTerms <- list()

lethalPhenotype <- c("preweaning lethality") %>% paste(., collapse = "|")
lethalRoot <- termDefinitions[grep(lethalPhenotype, termDefinitions$name), 'id']
lethalTerms <- get_descendants(ontologyFile, roots = paste(lethalRoot), exclude_roots = FALSE)


#brainPhenotype <- c("abnormal nervous system development") %>% paste(., collapse = "|")
brainPhenotype <- c("abnormal nervous system development") %>% paste(., collapse = "|")
brainRoot <- termDefinitions[grep(brainPhenotype, termDefinitions$name), 'id']
brainTerms <- get_descendants(ontologyFile, roots = paste(brainRoot), exclude_roots = FALSE)


kidneyPhenotype <- c("abnormal urinary system development") %>% paste(., collapse = "|")
kidneyRoot <- termDefinitions[grep(kidneyPhenotype, termDefinitions$name), 'id']
kidneyTerms <- get_descendants(ontologyFile, roots = paste(kidneyRoot), exclude_roots = FALSE)


liverPhenotype <- c("abnormal hepatobiliary system development") %>% paste(., collapse = "|")
liverRoot <- termDefinitions[grep(liverPhenotype, termDefinitions$name), 'id']
liverTerms <- get_descendants(ontologyFile, roots = paste(liverRoot), exclude_roots = FALSE)


listOfTerms[['Brain']] <- brainTerms
listOfTerms[['Kidney']] <- kidneyTerms
listOfTerms[['Liver']] <- liverTerms
listOfTerms[['embryonicLethal']] <- lethalTerms

devTissues <- c("Hindbrain", "Kidney", "Liver")


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

enrichmentList <- list()
termGenes <- list()
for (devTissue in devTissues) {
	load(glue("plsrModelForMedianPosEmbSplicingIn{devTissue}.rda"))
	testGeneList <- combinedLoadings[combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0.0001, ] %>% rownames() #%>% GeneID2entrez(.)
	backgroundList <- combinedLoadings[combinedLoadings$coefficients < 0.0, ] %>% rownames() #%>% GeneID2entrez(.)

	####select again top 100 genes####
	if (devTissue == "Kidney") {
		testGeneList <- combinedLoadings %>% .[.$FDR <= 0.05, ] %>% .[order(.$coefficients, decreasing = T), ] %>% rownames %>% .[1:45] #%>% GeneID2entrez(.)
	} else {
		#testGeneList <- combinedLoadings %>% .[.$FDR <= 0.05, ] %>% .[order(.$coefficients, decreasing = T), ] %>% rownames %>% .[1:100] %>% GeneID2entrez(.)
	}
	
	testGeneList <- testGeneList[testGeneList %in% allGenes]
	backgroundList <- backgroundList[backgroundList %in% allGenes]
	
	
	tissue <- gsub("Hindb", "B", devTissue)
	
	tempList <- list()
	for (tissue in names(listOfTerms)) {
		termsToUse <- listOfTerms[[tissue]] %>% unique
		geneList <- termMappings %>% .[.$MPID %in% paste(termsToUse), 'GeneSym'] %>% unique()
		
		csfInPhenotype <- testGeneList[testGeneList %in% geneList] 
		csfFraction <- length(csfInPhenotype) / length(testGeneList)
		ncsfInPhenotype <- backgroundList[backgroundList %in% geneList] 
		ncsfFraction <- length(ncsfInPhenotype) / length(backgroundList)
		
		tempDf <- fisherTestFunction(targets = geneList, testGeneList = testGeneList, backgroundList = backgroundList)
		tempList[[tissue]] <- tempDf
	}
	enrichmentList[[devTissue]] <- tempList
}

enrichmentDf <- lapply(enrichmentList, bind_rows, .id = "defectIn") %>% bind_rows(., .id = "devTissue")

enrichmentDf$devTissue <- gsub("Hindb", "B", enrichmentDf$devTissue)
colnames(enrichmentDf) <- c("factorsOf", "defectIn", "CSFoverlap", "nCSFoverlap", "oddsRatio", "pvalue")
rownames(enrichmentDf) <- NULL

enrichmentDf$FDR <- p.adjust(enrichmentDf$pvalue, method = "fdr")
enrichmentDf$nCSFoverlap <- -enrichmentDf$nCSFoverlap

plotDf <- enrichmentDf[enrichmentDf$pvalue < 0.1, ]
plotDf <- reshape2::melt(plotDf, id.vars = c("factorsOf", "defectIn", "oddsRatio", "pvalue", "FDR"))


ggTheme <- theme(axis.text.x = element_text(size = 14, angle = 0, vjust = 1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 14),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	#legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 14), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
	
plotDf$defectIn[c(1,5)] <- "Abnormal neuronal development"
plotDf$defectIn[c(2,3,4,6,7,8)] <- "preweaning lethality"
plotDf$variable <- gsub("overlap", "", plotDf$variable)
ggplot(data = plotDf, aes(x = defectIn, y = value, fill = variable)) + theme_bw() + ggTheme + 
    facet_wrap(~ factorsOf, scales = "free", ncol = 1, strip.position="right") +  geom_col() + coord_flip() 

ggsave(file = glue("mammPheaResultsUsingCuratedGeneSets.svg"), width =  8, height = 5) 
	
	

