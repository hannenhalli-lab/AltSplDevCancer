library(dplyr); library(magrittr); library(glue)
library(reshape2); library(ggplot2); library(data.table)
library(readxl);
library(clusterProfiler)
library(org.Hs.eg.db)
orgDb <- "org.Hs.eg.db"
library(tidyverse)
library(gplots)
library(pheatmap)

load("MFEnrichmentAnalysisOfBroadlyActiveDomains.Rda")
molecularEnrichmentList <- tissueEnrichmentList
names(molecularEnrichmentList) <- gsub("domainIntersect", "", names(molecularEnrichmentList)) %>% gsub("embPos", "EP-", .) %>% gsub("embNeg", "EN-", .)

load("BPEnrichmentAnalysisOfBroadlyActiveDomains.Rda")
biologicalEnrichmentList  <- tissueEnrichmentList
names(biologicalEnrichmentList) <- gsub("domainIntersect", "", names(biologicalEnrichmentList)) %>% gsub("embPos", "EP-", .) %>% gsub("embNeg", "EN-", .)


####following function used from https://stackoverflow.com/questions/30406560/multiple-intersection-of-lists####
listIntersectionFunction <- function(inputList1, inputList2) {
	#sapply(inputList1, function(x) sapply(inputList2, function(y) length(intersect(x,y))))
	sapply(inputList1, function(x) {
			normFactor1 = length(x); 
			sapply(inputList2, function(y) {
				normFactor2 = length(y); 
				normFactor <- normFactor1 * normFactor2;
				length(intersect(x,y)) / normFactor
				})
		})
}

geneSetOverlaps <- function(bpEnrichment, mfEnrichment) {
	biologicalGeneSets <- bpEnrichment %>% data.frame %>% subset(., select = c(Description, geneID)) %>% 
    				mutate(name=strsplit(geneID, "/")) %>% unnest(name) %>% subset(., select = c(Description,name)) %>% 
					split(., .$Description) %>% lapply(., function(x) x$name)
					
	molecularGeneSets <- mfEnrichment %>% data.frame %>% subset(., select = c(Description, geneID)) %>% 
					mutate(name=strsplit(geneID, "/")) %>% unnest(name) %>% subset(., select = c(Description,name)) %>% 
					split(., .$Description) %>% lapply(., function(x) x$name)

	overlapMatrix <- listIntersectionFunction(inputList1 = biologicalGeneSets, inputList2 = molecularGeneSets)
	return(overlapMatrix)
}

abbreviateFunction <- function(wordString, minwordLength = 5, cutoffLength = 6) {
	wordString <- strsplit(wordString, " ") %>% unlist %>% as.list()
	abbreviatedString <- lapply(wordString, function(x) {
		#y <- length(x)
		if (nchar(x) > cutoffLength) {
			abbreviate(x, minlength = minwordLength, method = "left.kept")
		} else {
			x
		}
	})
	unlist(abbreviatedString) %>% paste(., sep=" ", collapse = " ") %>% gsub(" +", " ", .)
}



overlapMatrixList <- list()
for (exonType in names(biologicalEnrichmentList)) {
	biologicalEnrichment <- biologicalEnrichmentList[[exonType]]
	molecularEnrichment <- molecularEnrichmentList[[exonType]]
	for (domainName in names(biologicalEnrichment)) {
		fileName = glue("functionalAssesmentFor{domainName}In{exonType}")
		bpEnrichment <- biologicalEnrichment[[domainName]]$enrichment %>% clusterProfiler::simplify(., cutoff = 0.7)
		mfEnrichment <- molecularEnrichment[[domainName]]$enrichment %>% clusterProfiler::simplify(., cutoff = 0.7)
		if (nrow(bpEnrichment) > 50) {
			bpEnrichment <- biologicalEnrichment[[domainName]]$level4 %>% clusterProfiler::simplify(., cutoff = 0.7)
		}
		#bpEnrichment@result$Description <- stringr::str_trunc(bpEnrichment@result$Description, 40)
		bpEnrichment@result$Description <- lapply(bpEnrichment@result$Description, abbreviateFunction) %>% unlist
		#bpEnrichment@result$Description <- abbreviate(bpEnrichment@result$Description, 1) %>% toupper
		
		#mfEnrichment@result$Description <- stringr::str_trunc(mfEnrichment@result$Description, 40)
		mfEnrichment@result$Description <- lapply(mfEnrichment@result$Description, abbreviateFunction) %>% unlist
		#mfEnrichment@result$Description <- abbreviate(mfEnrichment@result$Description, 1) %>% toupper
		
		if (nrow(bpEnrichment) == 0 & nrow(mfEnrichment) == 0) {
			overlapMatrix = NA
		}
		if (nrow(bpEnrichment) > 0 | nrow(mfEnrichment) > 0) {
			overlapMatrix = NA
		}
		if (nrow(bpEnrichment) > 0 & nrow(mfEnrichment) > 0) {
			overlapMatrix <- geneSetOverlaps(bpEnrichment = bpEnrichment, mfEnrichment = mfEnrichment)
		}
		overlapMatrixList[[fileName]] <- overlapMatrix
	}
}


lwid <- c(0.2, 0.9)
lhei <- c(0.3, 3.0, 0.6)
lmat = rbind(c(,3),c(2,1),c(0,4))


#library("colorspace")
#my_palette <- choose_palette()

my_palette <- colorRampPalette(c("khaki", "darkred"))(n = 10)
lwid <- c(0.25, 1.5)
lhei <- c(0.3, 3.0)
lmat = rbind(c(4,3),c(2,1))



for (fileName in names(overlapMatrixList)) {
	plotData <- overlapMatrixList[[fileName]]
	fileName <- glue("{fileName}-heatmap2")
	if (!is.na(plotData)) {
		png(filename = glue("{fileName}.tiff"), width = 20, height = 18, units = "in", pointsize = 14, bg = "white", res = 300, type = "quartz", antialias = 'none')
		rownames(plotData) <- stringr::str_trunc(rownames(plotData), 32)
		colnames(plotData) <- stringr::str_trunc(colnames(plotData), 32)
		#mfEnrichment@result$Description <- stringr::str_trunc(mfEnrichment@result$Description, 40)
		heatmap.2(plotData, trace = 'none', col = my_palette, lmat = lmat, lhei = lhei, lwid = lwid, key = T,  keysize = 0.5, key.title = NA, key.ylab = NA,
			margins = c(32,32), cexRow = 3, cexCol = 3, srtRow = 0, srtCol = 90, colsep=1:ncol(plotData), rowsep=1:nrow(plotData), sepcolor = "grey", sepwidth=c(0.01,0.01))
		#pheatmap(plotData, treeheight_row = 1, treeheight_col = 1, fontsize = 14, scale = 'row', angle_col = 90, color = my_palette)
		dev.off()
		if (ncol(plotData) > 1.5 * nrow(plotData)) {
			png(filename = glue("{fileName}.tiff"), width = 30, height = 18, units = "in", pointsize = 14, bg = "white", res = 300, type = "quartz", antialias = 'none')
			heatmap.2(plotData, trace = 'none', col = my_palette, lmat = lmat, lhei = lhei, lwid = lwid, key = F,  keysize = 0.5, key.title = NA,  key.ylab = NA,
				margins = c(32,32), cexRow = 3, cexCol = 3, srtRow = 0, srtCol = 90, colsep=1:ncol(plotData), rowsep=1:nrow(plotData), sepcolor = "grey", sepwidth=c(0.01,0.01))
			#pheatmap(plotData, treeheight_row = 1, treeheight_col = 1, fontsize = 14, scale = 'row', angle_col = 90, color = my_palette)
			dev.off()
		}
	}
}


