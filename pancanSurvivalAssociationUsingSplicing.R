library(dplyr); library(magrittr); 
library(preprocessCore); library(reshape2);
library(clusterProfiler); library(glue)
library(TCGAbiolinks); library(ggplot2)
library(survival); library(OneR)
library(ggpubr)
#library(survcomp)
naFilter <- function(dataFrame, cutoff) {
	nsize <- floor(ncol(dataFrame) * cutoff)
	indexes <- apply(dataFrame, 1, function(x) sum(is.na(x))) < nsize
	dataFrame[indexes, ]
}

GetBins <- function(x, binStep)
    {
		x = ceiling(x / binStep)
		x[x == 0] = 1
		return(x)
    }
	
	
	
#####https://medium.com/human-in-a-machine-world/custom-trycatch-in-r-b2106ff914a9#####
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}


survivalFunction <- function(clinicalData, splicing, test, binStep, binned) {
	#dataToTest = cbind(clinicalData, InclusionLevel = unlist(splicing))
	InclusionLevel = unlist(splicing)
	dataToTest = merge(clinicalData, InclusionLevel, by.x = "samples", by.y = 0)
	names(dataToTest)[names(dataToTest) == "y"] <- "InclusionLevel"
	dataToTest$age = scale(dataToTest$age,center = T,scale = T)
	if (test == "logRank") {
		fit <- NA; sd <- NA; groupMedians <- c(NA, NA) 
		medianLevel <- mean(dataToTest$InclusionLevel, na.rm = T, na.action = na.pass)
		dataToTest$BinnedInclusion <- ifelse(dataToTest$InclusionLevel >= medianLevel, "HighInclusion", "LowInclusion")
		groupMedians <- dataToTest %>% group_by(BinnedInclusion) %>% summarise(median = median(InclusionLevel, na.rm = T, na.action = na.pass))
		groupMedians <- na.omit(groupMedians)
		groupMedians <- groupMedians[ ,'median'] %>% unlist() %>% set_names(unlist(groupMedians[ ,'BinnedInclusion']))
		res <- list(fit = fit, sd = sd, groupMedians = groupMedians)
		tryCatch({
			fit <- survfit(Surv(time, status) ~ BinnedInclusion, data = dataToTest)
			sd <- survdiff(Surv(time, status) ~ BinnedInclusion, data = dataToTest)
			res <- list(fit = fit, sd = sd, groupMedians = groupMedians)
	    },error=function(x){})
		return(res)
	}
	
	if (test == "cox") {
		HR <- NA; Pval <- NA;
		res <- list(HR = HR, Pval = Pval)
		dataToTest$BinnedInclusion <- GetBins(dataToTest$InclusionLevel, binStep = binStep)
		tryCatch({
			if (binned) {
				cox.base = coxph(Surv(time,status) ~ BinnedInclusion + age, data = dataToTest)
				HR = summary(cox.base)$coefficients['BinnedInclusion',2]
				Pval = summary(cox.base)$coefficients['BinnedInclusion',5]
			} else {
				cox.base = coxph(Surv(time,status) ~ InclusionLevel + age, data = dataToTest)
				HR = summary(cox.base)$coefficients['InclusionLevel',2]
				Pval = summary(cox.base)$coefficients['InclusionLevel',5]
			}
			
			res <- list(HR = HR, Pval = Pval)
		},error=function(x){})
		return(res)
	}
}

processLogrankOutput <- function (survResults) {
	lowMedian <- lapply(survResults, function(x) x[["groupMedians"]][2]) %>% unlist()
	highMedian <- lapply(survResults, function(x) x[["groupMedians"]][1]) %>% unlist()
	groupMediansData <- data.frame(lowMedian = lowMedian, highMedian = highMedian)
	medianDifference <- groupMediansData[ ,2] - groupMediansData[ ,1]
	indexes <- !is.na(medianDifference)
	groupMediansData <- groupMediansData[indexes, ]

	survResults <- survResults[indexes]
	
	SurvFitData <- lapply(survResults, function(x) x[["fit"]])
	survDiffData <- lapply(survResults, function(x) x[["sd"]])
	
	#save(SurvFitData, SurvFitData, file = 'LogrankResultsForSurvivalWithDrosophilaLikeSetup.Rda')
	Pvalue <- lapply(survDiffData, function(x) (1 - pchisq(x$chisq, length(x$n) - 1))) %>% unlist() 
	EffectSize <- lapply(survDiffData, function(x) {ggVec = (x$obs / x$exp); ggVec[2] / ggVec[1]}) %>% unlist()
	
	survivalResults <- data.frame(Exon = names(survResults), EffectSize, Pvalue) %>% cbind(., groupMediansData)
}

getEvents <- function (AssociationClusters, ExonType, Association) {
	AssociationClusters <- AssociationClusters[AssociationClusters$medianExpression >= 0, ]
	AssociationClusters <- AssociationClusters[AssociationClusters$ExonType == "Embryonic", ]
	AssociationClusters$ExonType = paste(AssociationClusters$ExonType, AssociationClusters$Association, sep = ".")
	ExonType <- paste(ExonType, Association, sep = ".")
	events <- AssociationClusters[AssociationClusters$ExonType == ExonType, 'Exon'] %>% as.character()
	return(events)
}


cancerTypes <- c("LIHC", "GBM", "LGG", "KIRP", "BLCA", "ESCA", "LUSC", "PRAD", "READ", "STAD", "UCS", "CESC", "SKCM", "THCA", "UVM", "CHOL", "PAAD")
pathToExpressionFile <- "~/DataAndAnalysis/DevelopmentAndCaner/Splicing/Cancer/Kallisto/"

load("commonEventsAndFactorsInLiverBrainKidney.Rda")

devTissues <- c("Hindbrain", "Liver", "Kidney")
embPosList <- list()
for(devTissue in devTissues) {
	load(glue("../Normal/Kallisto/KallistoPSIpathwayCorrelationfor{devTissue}.Rda"))
	embPos <- getEvents(AssociationClusters = AssociationClusters, ExonType = "Embryonic", Association = "Positive")
	embPosList[[devTissue]] <- embPos
}

specificEvents <- embPosList %>% unlist() %>% table %>% .[. == 1] %>% names

allClinicalData <- read.table("/Users/singha30/DataAndAnalysis/SPAGE/RB1-Project/AveragedSamples/AllAvialableClinicalDataMod.txt", sep = "\t", header = T)

analysesType <- "cox"
commonSurvList <- list()
specificSurvList <- list()

for (cancerType in cancerTypes) {
	
	clinicalData <- allClinicalData[allClinicalData$type == cancerType, ]
	
	cancerFile <- glue("{pathToExpressionFile}TCGA/PSIvalues/SuppaPSIvaluesUsingKallistoIn{cancerType}.psi")
	cancerPSI <- read.table(paste(cancerFile), sep = "\t", header = T, check.names = F, na.strings = c(NA, "nan"))
	names(cancerPSI) <- gsub("\\.", "-", names(cancerPSI))

	cancerPSI <- naFilter(dataFrame = cancerPSI, cutoff = 0.50)

	cancerSamples <- TCGAquery_SampleTypes(names(cancerPSI), typesample = c("TP"))
	cancerPSI <- cancerPSI[, names(cancerPSI) %in% cancerSamples] %>% .[ , order(names(.))]

	names(cancerPSI) <- gsub("-01$", "", names(cancerPSI))
	clinicalData <- clinicalData[as.character(clinicalData$samples) %in% names(cancerPSI), ]
	clinicalData <- clinicalData[order(clinicalData$samples), ]
	cancerPSI <- cancerPSI[ ,names(cancerPSI) %in% clinicalData$samples]
	
	if(analysesType == "logRank") {
		
		InclusionLevel <- cancerPSI %>% .[rownames(.) %in% commonEvents, ] %>% apply(., 2, median, na.rm = T)
		survResultsList[[cancerType]] <- survivalFunction(clinicalData = clinicalData, splicing = InclusionLevel, test = analysesType, binStep = 10, binned = F)
	}else {
		if (doThis) {
			survResults <- list()
			for (exon in commonEvents) {
				exonIndex <- which(rownames(cancerPSI) == exon)
				#HR <- NA; Pval <- NA;
				res <- NA
				if (length(exonIndex) == 1) {
					InclusionLevel <- cancerPSI[exonIndex, ] * 100
					#res <- survivalFunction(clinicalData = clinicalData, splicing = InclusionLevel, test = "logRank", binStep = 10)
					res <- survivalFunction(clinicalData = clinicalData, splicing = InclusionLevel, test = "cox", binStep = 5, binned = T)
				}
				survResults[[exon]] <- res
				#HRdf[[i]] <- HR
				#pvalDf[[i]] <- Pval
				#if (i == 100) {
					#break
					#}
			}
			names(survResults) <- commonEvents
			EffectSize <- lapply(survResults, function(x) x[1]) %>% unlist()
			Pvalue <- lapply(survResults, function(x) x[2]) %>% unlist()
			survivalResults <- data.frame(Exon = names(survResults), EffectSize, Pvalue) %>% set_rownames(NULL)
			commonSurvList[[cancerType]] <- survivalResults
		}
		
		
		survResults <- list()
		for (exon in specificEvents) {
			exonIndex <- which(rownames(cancerPSI) == exon)
			#HR <- NA; Pval <- NA;
			res <- NA
			if (length(exonIndex) == 1) {
				InclusionLevel <- cancerPSI[exonIndex, ] * 100
				#res <- survivalFunction(clinicalData = clinicalData, splicing = InclusionLevel, test = "logRank", binStep = 10)
				res <- survivalFunction(clinicalData = clinicalData, splicing = InclusionLevel, test = "cox", binStep = 10, binned = T)
			}
			survResults[[exon]] <- res
			#HRdf[[i]] <- HR
			#pvalDf[[i]] <- Pval
			#if (i == 100) {
				#break
				#}
		}
		names(survResults) <- commonEvents
		EffectSize <- lapply(survResults, function(x) x[1]) %>% unlist()
		Pvalue <- lapply(survResults, function(x) x[2]) %>% unlist()
		survivalResults <- data.frame(Exon = names(survResults), EffectSize, Pvalue) %>% set_rownames(NULL)
		specificSurvList[[cancerType]] <- survivalResults
	}
	
}

save(commonSurvList, specificSurvList, file = "specificAndCommonEventsPancCanceSurvival.Rda")

commonResultsDf <- dplyr::bind_rows(commonSurvList, .id = "variable") %>% set_names(c("Cancer", "Exon", "EffectSize", "Pvalue")) %>% na.omit() %>% 
		split(., .$Cancer) %>% lapply(., function(x) {x %>% mutate(., fdr = p.adjust(.$Pvalue, method = "fdr"))}) %>% do.call("rbind", .) %>% 
		mutate(direction = ifelse(.$EffectSize > 1, "Worst survival", "Better survival")) #%>% subset(., fdr < 0.2) %>%

specificResultsDf <- dplyr::bind_rows(specificSurvList, .id = "variable") %>% set_names(c("Cancer", "Exon", "EffectSize", "Pvalue")) %>% na.omit() %>% 
		split(., .$Cancer) %>% lapply(., function(x) {x %>% mutate(., fdr = p.adjust(.$Pvalue, method = "fdr"))}) %>% do.call("rbind", .) %>%
		mutate(direction = ifelse(.$EffectSize > 1, "Worst survival", "Better survival")) #%>% subset(., fdr < 0.2) %>% 

commonTotalExons <- table(commonResultsDf$Cancer) %>% data.frame()
specificTotalExons <- table(specificResultsDf$Cancer) %>% data.frame()

commonResultsDf <- commonResultsDf %>% subset(., fdr < 0.30) %>% group_by(Cancer, direction) %>% summarize(., length(Exon))
specificResultsDf <- specificResultsDf %>% subset(., fdr < 0.30) %>% group_by(Cancer, direction) %>% summarize(., length(Exon))

cancerTypes <- commonResultsDf$Cancer %>% unique()
HazardType <- commonResultsDf$direction %>% unique()
tempDf <- expand.grid(cancerTypes, HazardType) %>% set_names(c("Cancer", "direction"))

commonResultsDf <- merge(tempDf, commonResultsDf, by = c("Cancer", "direction"), all.x = T) %>% set_names(c("cancer", "direction", "count")) %>% mutate(count = ifelse(is.na(.$count), 0, .$count))

cancerTypes <- specificResultsDf$Cancer %>% unique()
HazardType <- specificResultsDf$direction %>% unique()
tempDf <- expand.grid(cancerTypes, HazardType) %>% set_names(c("Cancer", "direction"))

specificResultsDf <- merge(tempDf, specificResultsDf, by = c("Cancer", "direction"), all.x = T) %>% set_names(c("cancer", "direction", "count")) %>% mutate(count = ifelse(is.na(.$count), 0, .$count))


commonResultsDf <- merge(commonResultsDf, commonTotalExons, by.x = "cancer", by.y = "Var1") %>% mutate(ExonType = "common EP events")
specificResultsDf <- merge(specificResultsDf, specificTotalExons, by.x = "cancer", by.y = "Var1") %>% mutate(ExonType = "specific EP events")

commonResultsDf$fraction <- commonResultsDf$count / commonResultsDf$Freq 
specificResultsDf$fraction <- specificResultsDf$count / specificResultsDf$Freq 

combinedDf <- rbind(commonResultsDf, specificResultsDf)

#commonResultsDf %>% mutate(count = ifelse(.$direction == "MoreHazard", 1, -1) * .$count) %>% group_by(cancer) %>% summarize(difference = sum(count))		 
#specificResultsDf %>% mutate(count = ifelse(.$direction == "MoreHazard", 1, -1) * .$count) %>% group_by(cancer) %>% summarize(difference = sum(count))		  

ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = 'none',
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
	
#dataForPlot <- combinedDf[combinedDf$ExonType == "common", ]
dataForPlot <- combinedDf
my_comparisons <- list(c("Better survival", "Worst survival"))
ggp <- ggplot(data= dataForPlot, aes(x = direction, y = fraction, col = direction)) + geom_boxplot() + facet_wrap(~ExonType) +
			scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) + theme_bw() + ggTheme + 
			stat_compare_means(method = "wilcox", comparisons = my_comparisons, size = 5, label.y = 0.5, tip.length = 0.1)
ggp + scale_color_manual(values=c("turquoise4", "violetred3")) 

fileName <- glue("survivalPlotForCommonEvents.svg")
ggsave(file = fileName, height = 3, width = 3.4)


dataForPlot <- combinedDf[combinedDf$direction == "Worst survival", ]
my_comparisons <- list(c("common", "specific"))
ggplot(data= dataForPlot, aes(x = ExonType, y = fraction)) + geom_boxplot() + facet_wrap(~direction) +
		scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) + theme_bw() + ggTheme + 
		stat_compare_means(method = "wilcox", comparisons = my_comparisons, size = 5, label.y = 0.5, tip.length = 0.1, paired = T)

fileName <- glue("survivalPlotForCommonEvents.svg")
ggsave(file = fileName, height = 4.7, width = 4.7)




#barcodes <- clinicalData$samples





