
getUniverse <- function(normalTissue, normalFilePath, allowedNA, numberOfSDs, testType, effectSize) {
	universeListIncrease <- list()
	universeListDecrease <- list()
	for (fileIndex in 1:length(normalTissue)) {
		tissueType <- normalTissue[fileIndex] %>% as.character()
		#tissueType <- normalTissue
		normalFile <- glue("{normalFilePath}{tissueType}.psi")
		normalPSI <- read.table(paste(normalFile), sep = "\t", header = T, check.names = F, na.strings = c(NA, "nan"))
		names(normalPSI) <- gsub("\\.", "-", names(normalPSI))
		
		normalNAcutoff <- (ncol(normalPSI) * (allowedNA / 100)) %>% floor()
		indexes <- apply(normalPSI, 1, function(x) sum(is.na(x))) < normalNAcutoff
		normalPSI <- normalPSI[indexes, ]
		
		
		if (testType == "outlier") {
			normalMean <- apply(normalPSI, 1, mean, na.rm = T)
			normalDeviation <- apply(normalPSI, 1, sd, na.rm = T)
		
			upperLimit <- normalMean + (numberOfSDs * normalDeviation)
			lowerLimit <- normalMean - (numberOfSDs * normalDeviation)
			
		
			#increasedUniverse <- table(upperLimit <= 1)['TRUE']
			#decreasedUniverse <- table(lowerLimit >= 0)['TRUE']
		
			#tempDf1 <- c(tissue = tissueType, eventType = "Increase", universe = unname(increasedUniverse))
			#tempDf2 <- c(tissue = tissueType, eventType = "Decrease", universe = unname(decreasedUniverse))
			
			increasedUniverse <- names(upperLimit[upperLimit < 1])
			decreasedUniverse <- names(lowerLimit[lowerLimit > 0])
			
			universeListIncrease[[fileIndex]] <- increasedUniverse
			universeListDecrease[[fileIndex]] <- decreasedUniverse
		}
		
		if (testType == "wilcox") {
			normalMean <- apply(normalPSI, 1, median, na.rm = T, na.action = na.pass)
			
			#increasedUniverse <- table(normalMean <= (1 - effectSize))['TRUE']
			#decreasedUniverse <- table(normalMean >= effectSize)['TRUE']
		
			#tempDf1 <- c(tissue = tissueType, eventType = "Increase", universe = unname(increasedUniverse))
			#tempDf2 <- c(tissue = tissueType, eventType = "Decrease", universe = unname(decreasedUniverse))
			
			increasedUniverse <- names(normalMean[normalMean <= (1 - effectSize)])
			decreasedUniverse <- names(normalMean[normalMean >= effectSize])
			
			universeListIncrease[[fileIndex]] <- increasedUniverse
			universeListDecrease[[fileIndex]] <- decreasedUniverse
			
		}
		}
		#tempDf <- rbind(tempDf1, tempDf2)
		#universeDf <- rbind(universeDf, tempDf)
    universeListIncrease <- unlist(universeListIncrease) %>% unique()
    universeListDecrease <- unlist(universeListDecrease) %>% unique()
		return(list(Increase = universeListIncrease, Decrease = universeListDecrease))
}

