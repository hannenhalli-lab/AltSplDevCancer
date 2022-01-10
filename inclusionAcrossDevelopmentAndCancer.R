library(dplyr); library(magrittr);
library(gtools); library(reshape2); library(ggplot2)
library(glue); 


cancerDevelopmentMapFunction <- function(devEvents, cancerEvents, cancerType, Tissue, eventType, cancerStages, expLevel, universeDf) {
	DifferentialCancerEvents <- cancerEvents %>% subset(., .$CancerType %in% cancerType & .$TissueType %in% Tissue)
	DifferentialCancerEvents <- DifferentialCancerEvents[DifferentialCancerEvents$CancerExpression >= expLevel, ]
	
	DifferentialCancerEvents <- DifferentialCancerEvents %>% subset(., select = -c(CancerType, TissueType, CancerExpression, NormalExpression, cancerMean, normalMean, normalDeviation)) %>% unique()
	
	exonByEventType <- split(DifferentialCancerEvents, DifferentialCancerEvents$CancerStage) %>% lapply(., function(x) {
					increasedExons = x[x$eventType == "Increase", 'Exon']; decreasedExons = x[x$eventType == "Decrease", 'Exon'];
					uniqueIncrease <- increasedExons[!(increasedExons %in% decreasedExons)];
					uniqueDecrease <- decreasedExons[!(decreasedExons %in% increasedExons)];
					df1 <- x[x$eventType == "Increase" & x$Exon %in% uniqueIncrease, ];
					df2 <- x[x$eventType == "Decrease" & x$Exon %in% uniqueDecrease, ];
					rbind(df1, df2) %>% return
					#return(list(increasedExons = increasedExons, decreasedExons = decreasedExons))
				})
	
	DifferentialCancerEvents <- do.call("rbind", exonByEventType)
	DifferentialCancerEvents <- DifferentialCancerEvents %>% .[.$eventType %in% eventType, ] %>% .[.$CancerStage %in% cancerStages, ]
	
	universe <- universeDf[universeDf$Universe == "Increase", 'Exon'] %>% as.character() %>% unique()
	events <- devEvents[devEvents$ExonType == "Embryonic.Positive", "Exon"]  %>% as.character() %>% unique()
	events1 <- events[events %in% universe]
	devEvents1 <- devEvents[devEvents$Exon %in% events1, ]
	
	universe <- universeDf[universeDf$Universe == "Decrease", 'Exon'] %>% as.character() %>% unique()
	events <- devEvents[devEvents$ExonType == "Embryonic.Negative", "Exon"]  %>% as.character() %>% unique()
	events2 <- events[events %in% universe]
	devEvent2 <- devEvents[devEvents$Exon %in% events2, ]
	
	devEvents <- rbind(devEvents1, devEvent2)
	
	backGroundInPathways <- table(devEvents$ExonType)
	
	DifferentialCancerEvents <- merge(DifferentialCancerEvents, devEvents[ ,c('Exon', 'ExonType')], by.x = "Exon", by.y = "Exon", all.x = T)
	indexes <- which(is.na(DifferentialCancerEvents$ExonType))
	DifferentialCancerEvents$ExonType[indexes] <- "Other"
	return(list(DifferentialCancerEvents = DifferentialCancerEvents, backGroundInPathways = backGroundInPathways, embPosInUniverse = devEvents1, embNegInUniverse = devEvent2))
}

DistributionPlotData <- function(InputFile, clusterFile) {
	InputFile$Exon <- rownames(InputFile)
	InputFile <- melt(InputFile)
	InputFile$ExonType = "Other"
	indexes <- match(InputFile$Exon, clusterFile$Exon)
	InputFile$ExonType <- clusterFile$ExonType[indexes]
	InputFile$ExonType  <- as.character(InputFile$ExonType)
	indexes <- is.na(InputFile$ExonType)
	InputFile$ExonType[indexes] <- "Other"
	return(InputFile)
}


getRandomExons <- function(normalFilePath, normalTissue, posExons, negExons, allowedNA, universeDf) {
	randomPosExons <- list()
	randomNegExons <- list()
	
	posSize <- length(posExons); negSize <- length(negExons);
	normalFile <- glue("{normalFilePath}{normalTissue}.psi")
	normalPSI <- read.table(paste(normalFile), sep = "\t", header = T, check.names = F, na.strings = c(NA, "nan"))
	names(normalPSI) <- gsub("\\.", "-", names(normalPSI))
	
	normalNAcutoff <- (ncol(normalPSI) * (allowedNA / 100)) %>% floor()
	indexes <- apply(normalPSI, 1, function(x) sum(is.na(x))) < normalNAcutoff
	normalPSI <- normalPSI[indexes, ]
	
	medianPSI <- normalPSI %>% apply(., 1, mean, na.rm = T, na.action = na.pass) ### take mean or median, either is fine
	
	increaseUniverse <- universeDf[universeDf$Universe == "Increase", 'Exon']
	decreaseUniverse <- universeDf[universeDf$Universe == "Decrease", 'Exon']
	
	posData <- normalPSI[medianPSI < 0.2, ] %>% .[rownames(.) %in% increaseUniverse, ] %>% rownames(.) %>% .[!(. %in% posExons)]
	negData <- normalPSI[medianPSI > 0.8, ]  %>% .[rownames(.) %in% decreaseUniverse, ] %>% rownames(.) %>% .[!(. %in% negExons)]
	
	set.seed(123)
	for(index in 1:100){
		posRan <- sample(posData, posSize)
		negRan <- sample(negData, negSize)
		
		randomPosExons[[index]] <- posRan
		randomNegExons[[index]] <- negRan
		
	}
	return(list(randomPositive = randomPosExons, randomNegtive = randomNegExons))
}

getEvents <- function (AssociationClusters, ExonType, Association) {
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	ExonType <- paste(ExonType, Association, sep = ".")
	events <- AssociationClusters[AssociationClusters$ExonType == ExonType, 'Exon'] %>% as.character()
	return(events)
}


fisherTestFunction <- function(list1, list2, universe) {
	list1 <- list1[list1 %in% universe]
	list2 <- list2[list2 %in% universe]
	
	universe <- length(universe)
	c1 <- sum(list1 %in% list2)
	c2 <- length(list1) - c1
	c3 <- length(list2) - c1
	c4 <- universe - (c3 + length(list1))
	
	fisherMatrix <- matrix(c(c1,c2,c3,c4), nrow = 2)
	fisherTest <- fisher.test(fisherMatrix)
	enrichment <- fisherTest$estimate
	pvalue <- fisherTest$p.value
	
	return(c(enrichment = enrichment, pvalue = pvalue))
}

#normalTissue = args[1]
#cancerType = args[2]
#devTissue = args[3]
#pseudoCount = args[4]

detectionMethod = "outlier"
if (detectionMethod == "outlier") {
	splicingEventsDf <- read.table("../Cancer/Kallisto/DifferentialEventsInCancerWithTCGAstagingUsingoutlierTestFromKallistoTrancripts-Comprehensive-2SD.txt", sep = "\t", header = T)
}



usedTissues <- c("Hindbrain", "Liver", "Kidney")
tissueParameters <- list()
tissueParameters[['Hindbrain']] <- list(devTissue = "HindBrain", normalTissue = c("BrainCerebellum", "BrainCortex"), cancerType =  c("GBM", "LGG"), cancerStages = c("Overall"))
tissueParameters[['Liver']] <- list(devTissue = "Liver", normalTissue = "Liver", cancerType =  c("LIHC"), cancerStages = c("Overall", "Early", "Late"))
tissueParameters[['Kidney']] <- list(devTissue = "Kidney", normalTissue = "KidneyCortex", cancerType =  c("KIRP", "KIRC"), cancerStages = c("Overall", "Early", "Late"))



source('getUniverse.R')
pathToExpressionFile <- "/Users/singha30/DataAndAnalysis/DevelopmentAndCaner/Splicing/Cancer/Kallisto/"
eventTypes <- c("Increase", "Decrease")
cancerStages <- c("Early", "Late", "Overall")
enrichmentValueList <- list()

numberOfSD <- 2
for (devTissue in usedTissues) {
	load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))
	
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType != "Not defined", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	
	#AssociationClusters <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Positive")
	
	
	#normalTissue <- "BrainCerebellum"
	
	cancerType <- tissueParameters[[devTissue]]['cancerType'] %>% unlist %>% unname
	normalTissue <- tissueParameters[[devTissue]]['normalTissue'] %>% unlist %>% unname
	
	
	normalFile <- glue("{pathToExpressionFile}Gtex/PSIvalues/SuppaPSIvaluesUsingKallistoIn")
	universeList <- getUniverse(normalTissue = normalTissue, normalFilePath = normalFile, allowedNA = 40, numberOfSDs = numberOfSD, testType = "outlier", effectSize = 0.2)
	universeDf <- stack(universeList) %>% set_colnames(c("Exon", "Universe"))
	
	DiffEventsOverlaps <- cancerDevelopmentMapFunction(devEvents = AssociationClusters, cancerEvents = splicingEventsDf, 
		cancerType = cancerType, Tissue = normalTissue, eventType = eventTypes, cancerStages = cancerStages, expLevel = 0, universe = universeDf)

	DifferentialCancerEvents <- DiffEventsOverlaps[['DifferentialCancerEvents']]
	backGroundInPathways <- DiffEventsOverlaps[['backGroundInPathways']]
	
	#universeDf <- data.frame(universeDf) %>% set_rownames(NULL)
	#universeDf$eventType <- as.character(universeDf$eventType)
	#universe <- ncol(psiPathwayCorrelation)
	enrichmentValues <- NULL
	for (stageIndex in 1:length(cancerStages)) {
		stage <- cancerStages[stageIndex]
		subData1 <- DifferentialCancerEvents[DifferentialCancerEvents$CancerStage == stage, ]
		
		if (nrow(subData1) > 100) {
			overlapStatics <- data.frame(table(subData1$eventType, subData1$ExonType)) %>% set_colnames(c("eventType", "ExonType", "Count"))
			backgroundInCancer <- table(subData1$eventType)
	
			for (eventIndex in 1:length(eventTypes)) {
				eventType <- eventTypes[eventIndex]
				subData2 <- overlapStatics[overlapStatics$eventType == eventType, ]
				universe <- universeDf[universeDf$Universe == eventType, 'Universe' ] %>% length()
				for(exonIndex in 1:length(names(backGroundInPathways))) {
					exonType <- names(backGroundInPathways)[exonIndex]
					c1 <- subData2[subData2$ExonType == exonType, 'Count']
					c2 <- backgroundInCancer[eventType] - c1
			
					#c3 <- sum(overlapStatics[overlapStatics$eventType != eventType & overlapStatics$ExonType == exonType, 'Count'])
					#c4 <- sum(overlapStatics[overlapStatics$eventType != eventType & overlapStatics$ExonType != exonType, 'Count'])
			
					c3 <- backGroundInPathways[exonType] - c1
					#c4 <- sum(backgroundInCancer[names(backgroundInCancer)!=eventType]) + sum(backGroundInPathways[names(backGroundInPathways)!=exonType])
					c4 <- universe - (c3 + backgroundInCancer[eventType])
			
					hyperP <- phyper(c1 - 1, backgroundInCancer[eventType], universe - backgroundInCancer[eventType], backGroundInPathways[exonType], lower.tail= FALSE)
			
					fisherMatrix <- matrix(c(c1,c2,c3,c4), nrow = 2)
					fisherTest <- fisher.test(fisherMatrix)
					enrichment <- fisherTest$estimate
					pvalue <- fisherTest$p.value
			
					standardError <- sqrt((1/c1) + (1/c2) + (1/c3) + (1/c4)) %>% format(.,digits = 3)
			
					if (exonType == "Embryonic.Positive") {
						exonType = "EP"
					}
					if (exonType == "Embryonic.Negative") {
						exonType = "EN"
					}
					fisherResult <- cbind(stage = stage, eventType = eventType, exonType = exonType, oddsRatio = unname(enrichment), pvalue = pvalue, standardError = unname(standardError), hyperGeometric = hyperP)
					enrichmentValues <- rbind(enrichmentValues, fisherResult)
				}
			}
		}
	}
	enrichmentValues <- data.frame(enrichmentValues) %>% set_rownames(NULL)
	enrichmentValues$fdr <- p.adjust(as.numeric(as.character(enrichmentValues$pvalue)), method = "fdr")
	enrichmentValues$oddsRatio <- as.numeric(as.character(enrichmentValues$oddsRatio))
	enrichmentValues$standardError <- as.numeric(as.character(enrichmentValues$standardError))
	levels(enrichmentValues$exonType) <- gsub("\\.", " ", levels(enrichmentValues$exonType))
	
	enrichmentValueList[[devTissue]] <- enrichmentValues
}
enrichmentValueDf <- bind_rows(enrichmentValueList, .id = "cancer")

#save(enrichmentValueDf, file = "enrichmentOfEmbryonicEventsInCancers.Rda")

stages <- c("Early", "Late")
dataForplot <- enrichmentValueDf[enrichmentValueDf$stage == "Overall", ]
dataForplot$cancer <- factor(dataForplot$cancer)
levels(dataForplot$cancer) <- gsub("Hindb", "B", levels(dataForplot$cancer))

barplot <- ggplot(data = dataForplot, aes(x = exonType, y = oddsRatio, color = exonType)) + theme_bw() + ylab("Odds Ratio") +
				geom_bar(stat="identity", fill="white") +  scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) + 
				geom_text(aes(label = formatC(fdr, format = "e", digits = 2),  size = 10, y = oddsRatio + 0.26), color = "black", position = position_dodge(0.9), size = 3.5)
								

ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
					
	
barplot + facet_grid(cancer ~ eventType, scale = "free_y") + ggTheme + scale_color_manual(values=c("turquoise4", "violetred3"))
ggsave(file = glue("plotForenrichmentOfEmbryonicEventsInCancers.svg"), width =  3.6, height = 5) 


stages <- c("Early", "Late")
dataForplot <- enrichmentValueDf[enrichmentValueDf$stage %in% stages, ]
dataForplot$cancer <- factor(dataForplot$cancer)
levels(dataForplot$cancer) <- gsub("Hindb", "B", levels(dataForplot$cancer))

barplot <- dataForplot %>% .[.$cancer == "Liver", ] %>% 
  ggplot(data = ., aes(x = exonType, y = oddsRatio, color = exonType)) + theme_bw() + ylab("Odds Ratio") +
  geom_bar(stat="identity", fill="white") +  scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) + 
  facet_grid(stage ~ eventType, scale = "free_y") + ggTheme + ggtitle("Liver cancer") + 
  geom_text(aes(label = formatC(fdr, format = "e", digits = 2),  size = 10, y = oddsRatio + 0.26), color = "black", position = position_dodge(0.9), size = 3.5)


ggsave(file = glue("plotForenrichmentOfEmbryonicEventsInLiver-EarlyLate.svg"), width =  3.6, height = 3.5) 

	
barplot <- dataForplot %>% .[.$cancer == "Kidney", ] %>% 
  ggplot(data = ., aes(x = exonType, y = oddsRatio, color = exonType)) + theme_bw() + ylab("Odds Ratio") +
  geom_bar(stat="identity", fill="white") +  scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) + 
  facet_grid(stage ~ eventType, scale = "free_y") + ggTheme +  ggtitle("Kidney cancer") + 
  geom_text(aes(label = formatC(fdr, format = "e", digits = 2),  size = 10, y = oddsRatio + 0.26), color = "black", position = position_dodge(0.9), size = 3.5)


ggsave(file = glue("plotForenrichmentOfEmbryonicEventsInKidney-EarlyLate.svg"), width =  3.6, height = 3.5) 




	
				##############FDR calculation using random embryonic events###############
				
				
posExons <- DiffEventsOverlaps[['embPosInUniverse']][ ,'Exon'] %>% as.character()
negExons <- DiffEventsOverlaps[['embNegInUniverse']][ ,'Exon'] %>% as.character()

randomEmbryonicExons <- getRandomExons(normalFilePath = normalFile[1], normalTissue = normalTissue[1], posExons = posExons, negExons = negExons, allowedNA = 40, universeDf = universeDf)


increasedUniverse <- universeDf %>% .[.$Universe == "Increase", 'Exon']; decreasedUniverse <- universeDf %>% .[.$Universe == "Decrease", 'Exon']

cancerIncreasedEvents <- splicingEventsDf %>% subset(., .$CancerType %in% cancerType & .$TissueType %in% normalTissue & .$eventType == "Increase", select = Exon) %>% 
						unlist() %>% as.character() %>% .[. %in% increasedUniverse]
						
cancerDecreasedEvents <- splicingEventsDf %>% subset(., .$CancerType %in% cancerType & .$TissueType %in% normalTissue & .$eventType == "Decrease", select = Exon) %>% 
						unlist() %>% as.character() %>% .[. %in% decreasedUniverse]


randomPositiveORs <- lapply(randomEmbryonicExons[['randomPositive']], function(x) fisherTestFunction(list1 = x, list2 = cancerIncreasedEvents, universe = increasedUniverse))
randomNegativeORs <- lapply(randomEmbryonicExons[['randomNegtive']], function(x) fisherTestFunction(list1 = x, list2 = cancerDecreasedEvents, universe = decreasedUniverse))

randomPositiveORs <- do.call("rbind", randomPositiveORs) %>% data.frame() %>% set_colnames(c("randomOddsratio", "pvalue")) %>% mutate(Exon = "EP")
randomNegativeORs <- do.call("rbind", randomNegativeORs) %>% data.frame() %>% set_colnames(c("randomOddsratio", "pvalue")) %>% mutate(Exon = "EN")

actualValue <- 3.5836571
ggplot(data = randomPositiveORs, aes(randomOddsratio)) + geom_density()

actualValue <- 3.6687271
ggplot(data = randomPositiveORs, aes(randomOddsratio)) + geom_density() #+ coord_cartesian(xlim = c(0.1, 4))


ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
                 axis.text.y = element_text(size = 12, face = "bold"),
                 strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
                 #legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
                 axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))



tempDf <- rbind(randomPositiveORs, randomNegativeORs) %>% data.frame()
ggplot(data = tempDf, aes(randomOddsratio, group = Exon, fill = Exon)) + geom_density() 




##### following is Obsolete for the current purpose###########
compareOddsratio = F
if (compareOddsratio) {
	oddsRatioComparisons <- NULL
	for (stageIndex in 1:length(cancerStages)) {
		stage <- cancerStages[stageIndex]
		subData1 <- enrichmentValues[enrichmentValues$stage == stage, ]
		for (exonIndex in 1:length(names(backGroundInPathways))) {
			exonType <- names(backGroundInPathways)[exonIndex]
			subData2 <- subData1[subData1$exonType == exonType, ]
		
			#####following steps taken from reference http://genometoolbox.blogspot.com/2014/06/test-for-difference-in-two-odds-ratios.html#####
			deltaOdds <- log(subData2[1, 'oddsRatio']) - log(subData2[2, 'oddsRatio'])
			standardError <- subData2[, 'standardError'] %>% .^2 %>% sum() %>% sqrt()
			zscore <- abs(deltaOdds) / standardError
			pvalue = 2*(1-pnorm(zscore))
		
			oddsComparison <- cbind(stage, exonType, event1 = as.character(subData2[1, 'eventType']), event2 = as.character(subData2[2, 'eventType']), deltaLogOdds = deltaOdds, standardError, zscore, pvalue)
			oddsRatioComparisons <- rbind(oddsRatioComparisons, oddsComparison)
		}
	}

	oddsRatioComparisons <- oddsRatioComparisons %>% data.frame %>% 
								mutate(fdr = p.adjust(as.numeric(as.character(.$pvalue)), method = "fdr")) %>%
								mutate(deltaLogOdds = as.numeric(as.character(.$deltaLogOdds))) %>%
								mutate(significance = factor(ifelse(.$fdr <= 0.1, "Yes","No"), levels = c("Yes", "No"))) %>%
								mutate(pvalue = as.numeric(as.character(.$pvalue)))

	oddsRatioComparisons$stage <- factor(oddsRatioComparisons$stage, levels = c("Overall", "Early", "Late"))
	ggDP <- ggplot(oddsRatioComparisons, aes(x=deltaLogOdds, y=reorder(exonType, deltaLogOdds), color = significance)) +
				   		geom_point(aes(size = -log10(pvalue))) + facet_wrap(~stage, scales = "fixed") +
						geom_segment(aes(x = 0, xend = deltaLogOdds, y = exonType, yend = exonType), linetype = 3, size = 1.4) +
						#scale_size_continuous(breaks = c(0.25, 0.75, 1.5, 2.0, 2.5)) + 
						scale_size_continuous(breaks = c(0.5, 2, 5, 10, 20)) + 
						ggtitle(glue("IncreaseVsDecrease: {cancerType}-{Tissue}")) +
						#scale_colour_gradient2(low = "red", mid = "white", high = muted("blue"), midpoint = 0.3) +
						xlab ("deltaLogOdds") +
						ylab ("exonType") + theme_bw ()
	
}		

clusteringBased = F
if(clusteringBased) {
	sortedPathways <- pathwayClusters[order(pathwayClusters$cluster), ] %>% rownames()
	colorNames <- c(group1 = "blue", group2 = "red")
	rowSidecolrs <- paste("group", pathwayClusters[order(pathwayClusters$cluster), 'cluster'], sep = "") %>% colorNames[.]

	NAsums <- apply(psiPathwayAssociation, 1, function(x) sum(!is.finite(x)))
	index <- which(NAsums <= 200)
	assocationMatrix <- psiPathwayAssociation[index, ] %>% 
							apply(., 1, function(x) {ifelse(is.finite(x), x, median(x))}) %>%
							.[match(sortedPathways, rownames(.)), ]


	###retain only protein coding genes####
	exonGenesDf <- 	exonGeneNames(exonList = colnames(assocationMatrix), GenesDf = GenesDf) %>% data.frame
	indexes <- which(exonGenesDf$geneType == "protein_coding")
	assocationMatrix <- assocationMatrix[ ,indexes]

	exonGenesDf <- 	exonGeneNames(exonList = colnames(assocationMatrix), GenesDf = GenesDf) %>% data.frame
	pcg <- exonGenesDf[exonGenesDf$geneType == "protein_coding", 'exonGenes'] %>% as.character()
	indexes <- match(pcg, rownames(normalGeneTPM))
	medianExpression <- apply(normalGeneTPM[indexes, ], 1, median, na.rm = T)

	indexes <- medianExpression > 1
	assocationMatrix <- assocationMatrix[ ,indexes]
	heatmap.2(as.matrix(assocationMatrix), trace = 'none', Rowv = F, scale = 'none', col = my_palette, cexRow = 0.7, cexCol = 0.7, RowSideColors = rowSidecolrs)
	legend("topright", legend = names(colorNames), col = colorNames, lty= 1, lwd = 5, cex=.8)


	distMat <- dist(t(assocationMatrix), method = "euclidean")
	exonClusters <- hclust(distMat, method = "complete" ) 
	#plot(exonClusters, cex = 0.6, hang = -1)
	clusters <- cutree(exonClusters, k = 2) 
	MappedexonClusters <- data.frame(colnames(assocationMatrix)) %>% mutate(cluster = clusters) %>% set_colnames(c("Exon", "Cluster"))
}

