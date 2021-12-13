#library(clusterProfiler)
library(dplyr); library(magrittr); library(ensembldb); library(EnsDb.Hsapiens.v86)
library(gtools); library(preprocessCore); library(reshape2); library(ggplot2)
library(clusterProfiler); library(msigdbr); library(gplots); library(FactoMineR); library(factoextra)
library(glue)
orgDb <- "org.Hs.eg.db"; EnsdB <- EnsDb.Hsapiens.v86

Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),

GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))

#############load transcript expression files and collpase the values by gene expression level#############


exonGeneNames <- function(exonList, GenesDf) {
	exonGenes <-  exonList %>% as.character() %>% gsub("\\..+", "", .)
	indexes <- match(exonGenes, GenesDf[ ,'GeneID'])
	exonGenes <- GenesDf[indexes, c('GeneName')] %>% as.character()
	geneType <- GenesDf[indexes, c('GeneType')] %>% as.character()
	return(cbind(exonGenes, geneType))
}

keggPathwayFile <- download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")
pathwayToGene <- keggPathwayFile[[1]] %>% set_colnames(c("PathwayID", "EntrezID"))
pathwayToName <- keggPathwayFile[[2]] %>% set_colnames(c("PathwayID", "PathwayName"))

hallmarks = F
if(hallmarks) {
	msigGenes = msigdbr(species = "Homo sapiens", category = "H")
	pathwayToGene = msigGenes %>% dplyr::select(gs_id, entrez_gene, gene_symbol) %>% as.data.frame() %>% set_colnames(c("PathwayID", "EntrezID", "GeneName"))
	pathwayToName = msigGenes %>% dplyr::select(gs_id, gs_name) %>% as.data.frame() %>% set_colnames(c("PathwayID", "PathwayName")) %>% unique()
}

convertedIds <- bitr(pathwayToGene$EntrezID, "ENTREZID", "ENSEMBL", orgDb, drop = F)
matches <- match(pathwayToGene$EntrezID, convertedIds$ENTREZID)
pathwayToGene$GeneIDs <- convertedIds[matches, 'ENSEMBL']

matches <- match(pathwayToGene$GeneIDs, GenesDf$GeneID)
pathwayToGene$GeneName <- GenesDf[matches, 'GeneName']

pathwaysToremove <- c("hsa01100","hsa01110","hsa01120","hsa01200","hsa01210","hsa01212","hsa01230","hsa01220") ###they are broad and redundant with others
pathwayToGene <- pathwayToGene[!(pathwayToGene$PathwayID %in% pathwaysToremove), ]
####load DNA sequence if needed######
#dna <- ensembldb:::getGenomeTwoBitFile(EnsdB)


pathwayAnnotationFunction <- function(cosineMatrix, embryonicOffset, K) {
	#######pathway separation####
	distMat <- dist(cosineMatrix, method = "euclidean")
	hrcCluster <- hclust(distMat, method = "complete" )
	
	mainClusters <- c("Embryonic", "Adult")
	clusterList <- list()
	clusters <- cutree(hrcCluster, k = K)
	MappedClusters <- data.frame(cosineMatrix) %>% mutate(cluster = clusters) %>% set_rownames(rownames(cosineMatrix))

	adultOffset <- ncol(cosineMatrix) - embryonicOffset
	stagevector <- c(rep("Embryonic", embryonicOffset), rep("Adult", adultOffset)) %>% set_names(colnames(cosineMatrix))

	meltedClusters <- melt(MappedClusters, id.vars = "cluster")
	meltedClusters$variable <- gsub("X", "", meltedClusters$variable)
	indexes <- match(meltedClusters$variable, names(stagevector))
	meltedClusters$stage <- stagevector[indexes]
	
	subData <- meltedClusters[meltedClusters$cluster == 1, ]
	
	medianValues <- meltedClusters %>% subset(., cluster == 1)  %>%
				group_by(stage) %>% 
				mutate(median = median(value, na.rm = T)) %>%
				dplyr::select(stage, median) %>% unique()
				
	index = which(medianValues$median == max(medianValues$median))
	cluster1 = medianValues$stage[index]
	clusterList[[1]] <- cluster1
	
	remainingCluster <- mainClusters[!(mainClusters %in% cluster1)]
	#index2 = which(medianValues$median == min(medianValues$median))
	
	medianValues <- meltedClusters %>% subset(., cluster != 1)  %>%
				group_by(cluster) %>% 
				mutate(median = median(value, na.rm = T)) %>%
				dplyr::select(cluster, median) %>% unique()
	
	index = which(medianValues$median == max(medianValues$median))
	cluster2 = medianValues$cluster[index]
	clusterList[[cluster2]] <- remainingCluster
	
	clusterAnnotations <- c(1,2) %>% set_names(unlist(clusterList))
	
	if (K == 3) {
		middleClusterNumber <- lapply(clusterList, is.null) %>% unlist %>% which()
		clusterList[[middleClusterNumber]] <- "Middle"
		#cluster2 = medianValues$stage[index2]
		clusterAnnotations <- c(1,2,3) %>% set_names(unlist(clusterList))
	}
	
	clusterNameVec <- names(clusterAnnotations[MappedClusters$cluster])
	clusterAnnotations <- data.frame(pathwayId = rownames(MappedClusters), cluster = MappedClusters$cluster, clusterType = clusterNameVec)
	#clusterAnnotations$clusterType = as.character(clusterAnnotations$clusterType)
	#indexes <- is.na(clusterAnnotations$clusterType)
	#clusterAnnotations$clusterType[indexes] <- "Middle"
	return(clusterAnnotations)
}


exonAnnotationFunction <- function(psiPathwayAssociation, psiPathwaySignificance, pathwayClusters, AssociationType, fdrLevel, K) {
	fdr <- apply(psiPathwaySignificance, 1, p.adjust, method = "fdr") %>% t() %>% data.frame
	significanceMatrix <- apply(fdr, 2, function(x) ifelse(x <= fdrLevel, 1, 0))
	
	psiPathwayAssociation <- psiPathwayAssociation * significanceMatrix
	psiPathwayAssociation$Exon <- rownames(psiPathwayAssociation)
	psiPathwayAssociation <- melt(psiPathwayAssociation)
	names(psiPathwayAssociation) <- c("Exon", "Pathway", "Association")
	
	
	psiPathwayAssociation$PathwayType <- "Other"
	matches <- match(psiPathwayAssociation$Pathway, pathwayClusters$pathwayId)
	psiPathwayAssociation$PathwayType <- pathwayClusters$clusterType[matches]
	clusterSizes <- table(pathwayClusters$clusterType)
	
	if (AssociationType == "Positive") {
		Associations <- psiPathwayAssociation[psiPathwayAssociation$Association > 0, ] %>% na.omit()
	}
	if (AssociationType == "Negative") {
		Associations <- psiPathwayAssociation[psiPathwayAssociation$Association < 0, ] %>% na.omit()
	}
	
	AssociationFrequency <- data.frame(table(Associations$Exon, Associations$PathwayType)) %>% set_colnames(c("Exon","Cluster", "Count"))
	Embryonic <- AssociationFrequency[AssociationFrequency$Cluster == "Embryonic", ] %>% .[order(.$Exon), ]
	Middle <- AssociationFrequency[AssociationFrequency$Cluster == "Middle", ] %>% .[order(.$Exon), ]
	Adult <- AssociationFrequency[AssociationFrequency$Cluster == "Adult", ] %>% .[order(.$Exon), ]

	Embryonic$FractionSignificant <- Embryonic$Count / clusterSizes['Embryonic']
	Middle$FractionSignificant <- Middle$Count / clusterSizes['Middle']
	Adult$FractionSignificant <- Adult$Count / clusterSizes['Adult']
	
	if (K == 3) {
		AssociationClusters <- data.frame(Exon = as.character(Embryonic$Exon), Embryonic = Embryonic$FractionSignificant, Middle = Middle$FractionSignificant, Adult = Adult$FractionSignificant, Association = AssociationType)
	}else{
		AssociationClusters <- data.frame(Exon = as.character(Embryonic$Exon), Embryonic = Embryonic$FractionSignificant, Adult = Adult$FractionSignificant, Association = AssociationType)
	}
	
	#ggplot(data = clustersForPlot, aes(x = cluster1, y = cluster2)) + theme_bw() + geom_point(color="grey") + facet_wrap(~Association)
	
	AssociationClusters$ExonType = "Not defined"
	indexes <- AssociationClusters$Embryonic >= 0.1 & AssociationClusters$Adult <= 0.05
	AssociationClusters$ExonType[indexes] <- "Embryonic"
	indexes <- AssociationClusters$Adult >= 0.1 & AssociationClusters$Embryonic <= 0.05
	AssociationClusters$ExonType[indexes] <- "Adult"
	
	if (K == 3) {
		indexes <- AssociationClusters$Middle >= 0.1 & AssociationClusters$Embryonic <= 0.05 & AssociationClusters$Adult <= 0.05
		AssociationClusters$ExonType[indexes] <- "Middle"
	}
	#backGroundInPathways <- table(AssociationClusters$ExonType)
	return(AssociationClusters)
}



Tissue <- "Kidney"; cancerType <-  "KIRP"

expressionFile <- glue("development/Expression/KallistoTPMvaluesforTranscriptExpressionIn{Tissue}")
psiFile <- glue("development/PSIvalues/SuppaPSIvaluesUsingKallistoIn{Tissue}.psi")



#####get expression data####
normalExpression <- read.table(expressionFile, sep = "\t", header = T, check.names = F)
rownames(normalExpression) <- gsub("\\.\\d+", "", rownames(normalExpression))
matches <- match(rownames(normalExpression), Transcripts$tx_id)
Genes <- Transcripts[matches, 'gene_name']
normalExpression <- cbind("GeneName" = Genes, normalExpression)
normalGeneTPM <- aggregate(.~GeneName, data = normalExpression, sum, na.rm = T, na.action = na.pass)

normalGeneTPM <- normalGeneTPM[ ,-1] %>% as.matrix() %>% normalize.quantiles() %>% 
					set_rownames(normalGeneTPM$GeneName) %>% 
					set_colnames(colnames(normalGeneTPM)[-1]) %>% t() %>% data.frame() %>%
					mutate(Stage = gsub("_Rep\\d+", "", colnames(normalGeneTPM)[-1])) %>% 
					aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
					set_rownames(.$Stage) %>% .[,-1] %>% t() %>% data.frame() %>%
					set_rownames(normalGeneTPM$GeneName)
					
allPathways <- unique(pathwayToGene$PathwayID)


medianPathwayExpression <- NULL
pathwayNames <- NULL
for (index in 1:length(allPathways)) {
	pathway <- allPathways[index]
	pathwayNames[index] <- pathway
	pathwayGenes <- pathwayToGene[pathwayToGene$PathwayID == pathway, 'GeneName'] %>% unique() %>% as.character()
	subData <- normalGeneTPM[rownames(normalGeneTPM) %in% pathwayGenes, ]
	medianExpression <- apply(subData, 2, median, na.rm = T, na.action = na.pass)
	#medianExpression <- apply(subData, 1, function(x) geneRank(x))
	medianPathwayExpression <- cbind(medianPathwayExpression,  medianExpression)
}
medianPathwayExpression <- data.frame(medianPathwayExpression) %>% set_colnames(pathwayNames)


####heatmap clustering#####
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 40)
indexes <- match(mixedsort(rownames(medianPathwayExpression)), rownames(medianPathwayExpression))
medianPathwayExpression <- medianPathwayExpression[indexes, ]

indexes <- which(colSums(medianPathwayExpression) != 0)
medianPathwayExpression <- medianPathwayExpression[ ,indexes]
heatmap.2(t(normalize.quantiles(as.matrix(medianPathwayExpression))), dendrogram = 'row', Colv = F, trace = 'none', scale = 'none', col = my_palette)


forPCA <-  medianPathwayExpression %>% as.matrix() %>% normalize.quantiles() %>% 
					set_rownames(rownames(medianPathwayExpression)) %>% 
					set_colnames(colnames(medianPathwayExpression)) %>% data.frame()
					
rownames(forPCA) <-  gsub(glue("{Tissue}_"), "", rownames(forPCA))
dropColumns <- !(apply(forPCA, 2, function(x) sum(is.na(x))) == nrow(medianPathwayExpression))
forPCA <- forPCA[ ,dropColumns]

pathwaysPCA <- PCA(forPCA, graph = FALSE, ncp = ncol(forPCA))
components <- data.frame(pathwaysPCA$ind$coord)

components$Stage <- gsub(glue("{Tissue}_|days"), "", rownames(components)) %>% as.character() %>% as.numeric()
#rownames(components) <-  gsub(glue("{Tissue}_"), "", rownames(components))
ggplot(components, aes(x=Dim.1, y=Dim.2)) + geom_point(size = 4) + scale_size_continuous(breaks = c(0.25, 0.75, 1.5, 2.0, 2.5))

#library(factoextra)
#fviz_pca(pathwaysPCA,  axes = c(1,2), label = "ind", invisible="none")

######Get loadings#######
componentsToUse <- 5
loadingsDf <- pathwaysPCA$var$coord
loadingsForCosine <- loadingsDf[ ,1:componentsToUse]
contibution <- abs(loadingsForCosine) %>% sweep(., 2, colSums(.), "/")
#loadingsForCosine <- loadingsForCosine[contibution[,1] > 0.002, ]
scoresForCosine <- components[ ,1:componentsToUse]
componentWeights <- pathwaysPCA$eig[,2][1:componentsToUse]

#library(proxy)
#cosineSimilarity <- simil(as.matrix(loadingsForCosine), as.matrix(scoresForCosine), method = "cosine")
#heatmap.2(t(cosineSimilarity), dendrogram = 'col', Rowv = F, trace = 'none', scale = 'none', col = my_palette)
#heatmap.2(cosineSimilarity, dendrogram = 'row', Rowv = T, Colv = F, trace = 'none', scale = 'none', col = my_palette)


###########################Weighted cosine similarity###########################
#####cosine similarity formula is modified from######
###https://stats.stackexchange.com/questions/31565/is-there-an-r-function-that-will-compute-the-cosine-dissimilarity-matrix
cosineSimilarity <- function(index, loadingData, scoreData, weights, contribution, useContribution) 
{
    scaledWeights <- weights / sum(weights)
	A = loadingData[index, ]
	cGG = contibution[index, ]
	cosineData <- NULL
	for (scoreIndex in 1:nrow(scoreData)) {
		 B = scoreData[scoreIndex, ]
		 if (useContribution) {
		 	cosineData[scoreIndex] <- sum(A * B * cGG * scaledWeights )/sqrt(sum(A^2)*sum(B^2))
		 }else{
		 	cosineData[scoreIndex] <- sum(A * B * scaledWeights )/sqrt(sum(A^2)*sum(B^2))
		 }
		 #cosineData[scoreIndex] <- sum(A * B )/sqrt(sum(A^2)*sum(B^2))
	}
    return(cosineData)
}

weightedCosine <- sapply(1:nrow(loadingsForCosine), cosineSimilarity, loadingData = loadingsForCosine, scoreData = scoresForCosine, weights = componentWeights, contribution = contribution, useContribution = F)
rownames(weightedCosine) <- rownames(scoresForCosine)
colnames(weightedCosine) <- rownames(loadingsForCosine)
weightedCosine <- weightedCosine %>% t()
heatmap.2(weightedCosine, dendrogram = 'row', Rowv = T, Colv = F, trace = 'none', scale = 'none', col = my_palette)


pathwayClusters <- pathwayAnnotationFunction(cosineMatrix = weightedCosine, embryonicOffset = 13, K = 2)

			
CorrelationFunction <- function(index, PSIfile, pathwayFile) {
	exonUsage <- PSIfile[index, ] %>% unlist()
	corList <- list()
	pvalList <- list()
	outputList <- list()
	for(pathwayIndex in 1:ncol(pathwayFile)) {
		pathwayUsage <- pathwayFile[ ,pathwayIndex]
		correlations <- cor.test(exonUsage, pathwayUsage)
		corList[pathwayIndex] <- correlations$estimate
		pvalList[pathwayIndex] <-  correlations$p.value
	}
	#outputList[['corr']] <- unlist(corList)
	#outputList[['pval']] <- unlist(pvalList)
	#return(outputList)
	cbind(correlation = unlist(corList), pvalue = unlist(pvalList))
}



#####PSI file and PSI-pathway correlation#####
devPSIfile <- read.table(psiFile, sep = "\t", header = T, na.strings = c('nan', "NA"), check.names = F)
#tpmcutOff <- 1
NAsums <- apply(devPSIfile,1, function(x) sum(is.na(x)))
NAthreshold <- ceiling(ncol(devPSIfile) / 2)
indexes <- which(NAsums < NAthreshold)
devPSIfile <- devPSIfile[indexes, ]

transposednormalFile <- data.frame(t(devPSIfile)) %>% 
							mutate(Stage = gsub("_Rep\\d+", "", names(devPSIfile))) %>% 
							aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
							set_rownames(.$Stage) %>% 
							.[,-1] %>% t() %>% data.frame() %>%
							set_rownames(rownames(devPSIfile))
			
transposednormalFile <- transposednormalFile[, mixedsort(names(transposednormalFile))]	



psiPathwayCorrelation <- sapply(1:nrow(transposednormalFile), CorrelationFunction, PSIfile = transposednormalFile, pathwayFile = medianPathwayExpression)


psiPathwayAssociation <- psiPathwayCorrelation[1:ncol(medianPathwayExpression), ] %>% t() %>% data.frame %>%
				set_rownames(rownames(transposednormalFile)) %>% set_colnames(colnames(medianPathwayExpression))
								
psiPathwaySignificance <- psiPathwayCorrelation[(ncol(medianPathwayExpression) + 1):nrow(psiPathwayCorrelation), ]  %>% t() %>% data.frame %>%
				set_rownames(rownames(transposednormalFile)) %>% set_colnames(colnames(medianPathwayExpression))					


####retain only protein coding genes####
exonGenesDf <- 	exonGeneNames(exonList = rownames(psiPathwayAssociation), GenesDf = GenesDf) %>% data.frame
indexes <- which(exonGenesDf$geneType == "protein_coding")
psiPathwayAssociation <- psiPathwayAssociation[indexes, ]
psiPathwaySignificance <- psiPathwaySignificance[indexes, ]


#tempDf <- significanceMatrix[indexes, ]

positiveAssociationCluster <- exonAnnotationFunction(psiPathwayAssociation = psiPathwayAssociation, psiPathwaySignificance = psiPathwaySignificance, 
						pathwayClusters = pathwayClusters, AssociationType = "Positive", fdrLevel = 0.05, K = 2)
	
negativeAssociationCluster <- exonAnnotationFunction(psiPathwayAssociation = psiPathwayAssociation, psiPathwaySignificance = psiPathwaySignificance, 
						pathwayClusters = pathwayClusters, AssociationType = "Negative", fdrLevel = 0.05, K = 2)	

AssociationClusters <- rbind(positiveAssociationCluster, negativeAssociationCluster)
exonGenesDf <- 	exonGeneNames(exonList = AssociationClusters$Exon, GenesDf = GenesDf) %>% data.frame
indexes <- match(exonGenesDf$exonGenes, rownames(normalGeneTPM))
medianExpression <- apply(normalGeneTPM[indexes, ], 1, median, na.rm = T)
AssociationClusters <- AssociationClusters %>% mutate(geneName = rownames(normalGeneTPM)[indexes], medianExpression = medianExpression)
					
save(medianPathwayExpression, pathwaysPCA, pathwayClusters, psiPathwayCorrelation, AssociationClusters, file = glue("KallistoPSIpathwayCorrelationfor{Tissue}.Rda"))
