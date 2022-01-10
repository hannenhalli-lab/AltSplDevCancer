#CNOT6,FAM98B,GTPBP2,HNRNPA1,HNRNPA3,HNRNPD,HNRNPM,KHSRP,NCL,PABPC1,PTBP1,RBM3,RBM38,SF3A1,SF3A2,SF3A3,SF3B3,SNRPB,SRPK1,TTF2,YBX1

#library(clusterProfiler)
library(dplyr); library(magrittr); library(gplots); library(glue)
library(gtools); library(reshape2); library(ggplot2); library(TFEA.ChIP)
library(ensembldb); library(EnsDb.Hsapiens.v86); EnsdB <- EnsDb.Hsapiens.v86
source("getEnrichedGenes.R")

Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),

GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))

#database = "RM_TSS5kb.Rdata"
database = "RM_TSS1kb.Rdata"
#database = "default"
if (database == "default") {
}else {
	load(glue("TFdatabases/{database}"))
	MetaData$Accession <- MetaData$Name
	commonNames <- colnames(Mat01)[colnames(Mat01) %in% MetaData$Accession]
	Mat01 <- Mat01[,colnames(Mat01) %in% commonNames]
	MetaData <- MetaData[MetaData$Accession %in% commonNames, ]
}

entrezToGeneName <- function(inputIDs, mapTable) {
	indexes <- match(inputIDs, mapTable$ENTREZ.ID)
	mapTable[indexes, 'GENE.ID']
}

##### following is obsolete#####
getTargetGenes <- function(testGeneList, backgroundList, pvalueDf) {
	#geneLists <- getGeneList(testGeneList,backgroundList) #generates list genes with binding in test set
	geneLists <- getGeneList(backgroundList, testGeneList) #generates list genes with binding in test set
	names(geneLists) <- get_chip_index()$TF
	tfNames <- data.frame(Accession = pvalueDf$Accession, TF = pvalueDf$TF)
	#indexes <- match(pvalueDf$Accession, names(geneLists))
	#geneLists <- geneLists[indexes]
	geneNames <- lapply(geneLists, function(inputIDs) entrezToGeneName(inputIDs, mapTable))
	geneNames <- split(unlist(geneNames, use.names = FALSE), rep(names(geneNames), lengths(geneNames))) %>%
				lapply(., function(x) unique(x))
}

tcgaStaging <- function(cancerType, samples) {
	stageInformation <- "Not defined"
	clinical <- GDCquery_clinic(paste("TCGA-",cancerType, sep = ""), type = "clinical")
	if (is.null(clinical$ajcc_pathologic_stage)) {
		stageInformation = "Not defined"
		if (cancerType == "OV") {
			clinical <- data.frame(Barcode = clinical$submitter_id, Stage = clinical$figo_stage)#, MetaStasis = clinical$ajcc_pathologic_m)
			clinical$Stage <- gsub("[A-HJ-UW-Z]$", "", clinical$Stage)
			if (length(table(clinical$Stage)) > 0) {
				stageInformation = "defined"
			}
		}
	} else {
		clinical <- data.frame(Barcode = clinical$submitter_id, Stage = clinical$ajcc_pathologic_stage)#, MetaStasis = clinical$ajcc_pathologic_m)
		clinical$Stage <- gsub("[A-HJ-UW-Z]$", "", clinical$Stage)
		if (length(table(clinical$Stage)) > 0) {
			stageInformation = "defined"
		}
	}
	
	samples <- gsub("-\\w+$", "", samples)
	clinical <- clinical[match(samples, clinical$Barcode), ]

	early = c("Stage I", "Stage II"); late = c("Stage III", "Stage IV")
	earlySamples <- clinical[clinical$Stage %in% early, 'Barcode'] %>% as.character()
	lateSamples <- clinical[clinical$Stage %in% late, 'Barcode'] %>% as.character()
	stageDf <- rbind(data.frame(sampleName = earlySamples, sampleType = "Early Cancer"), data.frame(sampleName = lateSamples, sampleType = "Late Cancer"))
	
	return(stageDf)
}


#####load all splicing factors
splicingFactors <- read.table("/Users/singha30/DataAndAnalysis/ExonicRepeatsInCancer/SplicingFactorsGO-0008380fromAmigo", sep = "\t", header = F)


devTissue <- "Liver"
load(glue('plsrModelForMedianPosEmbSplicingIn{devTissue}.rda'))
#####load results of correlation and regression analysis
#####load development and cancer expression files
pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/Normal/Kallisto/"

comparison <- "headTohead"
if (comparison == "headTohead") {
	testGeneList <- combinedLoadings[combinedLoadings$FDR <= 0.05 & combinedLoadings$coefficients > 0.0001, ] %>% rownames() %>% GeneID2entrez(.)
	backgroundList <- combinedLoadings[combinedLoadings$coefficients < 0.0, ] %>% rownames() %>% GeneID2entrez(.)
	
	cmMatrix <- contingency_matrix(testGeneList,backgroundList) #generates list of contingency tables, one per dataset
	pvalueDf <- getCMstats(cmMatrix) #generates list of p-values and OR from association test
	rankedTFs <- rankTFs(na.omit(pvalueDf), rankMethod = "wilcoxon", makePlot = F)
}

#save(pvalueDf, rankedTFs, combinedLoadings, file = glue("TFbindingEnrichmentIn1kbPromotersUsingTFEA.ChIPremapDataIn{devTissue}SplicingFactors.Rda"))
save(pvalueDf, rankedTFs, combinedLoadings, file = glue("TFbindingEnrichmentIn1kbPromotersUsingTFEA.ChIPremapDataIn{devTissue}SplicingFactors-NegativeAsBackground.Rda"))

