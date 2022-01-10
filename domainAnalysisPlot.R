library(dplyr); library(magrittr); library(glue)
library(reshape2); library(ggplot2); library(data.table)
library(readxl)    

####following function from https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames
read_excel_allsheets <- function(filename, tibble = FALSE) {
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

input = "textFiles"

if (input == "excel") {
	fileName <- "ExonEnrichedProteinDomains.xlsx"
	domainEnrichment <- read_excel_allsheets(fileName)
	significantDomainNames <- lapply(domainEnrichment, function(x) x[x$padj < 0.1, 'domain']) %>% unlist() %>% unique()
	significantDomainsOR <- lapply(domainEnrichment, function(x) x[x$domain %in% significantDomainNames, c('domain', 'OR')])
	significantDomainsFDR <- lapply(domainEnrichment, function(x) x[x$domain %in% significantDomainNames, c('domain', 'padj')])
	
}

if (input == "textFiles") {
	files <- "domainEnrichmentMergedExons/enrichment/"
	listOfFiles <- list.files(files)
	domainEnrichment <- lapply(listOfFiles, function(x) read.table(glue("{files}{x}"), sep = "\t", header = T))
	names(domainEnrichment) <- gsub("Vs.+", "", listOfFiles)
	significantDomainNames <- lapply(domainEnrichment, function(x) x[x$padj < 0.1 & x$groupFreq >= 3, 'domain']) %>% unlist() %>% unique()
	significantDomainsOR <- lapply(domainEnrichment, function(x) x[x$domain %in% significantDomainNames, c('domain', 'OR')])
	significantDomainsFDR <- lapply(domainEnrichment, function(x) x[x$domain %in% significantDomainNames, c('domain', 'padj')])
	
}

significantDomainsOR <- bind_rows(significantDomainsOR, .id = "exonType")
significantDomainsFDR <- bind_rows(significantDomainsFDR, .id = "exonType")

significantDomainsOR$FDR <- significantDomainsFDR$padj
significantDomainsOR$significance <- ifelse(significantDomainsFDR$padj < 0.1, "Significant", "Not significant")
significantDomainsOR$Tissue <- gsub("emb|Pos|Neg", "", significantDomainsOR$exonType)
significantDomainsOR$exonType <- gsub("Liver|Brain|Kidney", "", significantDomainsOR$exonType)
significantDomainsOR$enrichment <- ifelse(significantDomainsOR$OR > 1, "enriched", "depeleted")
significantDomainsOR$OR <- as.numeric(as.character(significantDomainsOR$OR))
significantDomainsOR$FDR <- as.numeric(as.character(significantDomainsOR$FDR))


significantDomainsOR[is.infinite(significantDomainsOR$OR), 'OR'] <- 150

#significantDomainsOR<- mutate(significantDomainsOR, xLines = as.numeric(factor(Tissue, levels = unique(Tissue))),  xLines = xLines + 0.5)
significantDomainsOR<- mutate(significantDomainsOR, xLines = 1.5)		
significantDomainsOR <- mutate(significantDomainsOR, yLines = as.numeric(factor(domain, levels = unique(domain))),  yLines = yLines + 0.5)
significantDomainsOR$exonType <- factor(significantDomainsOR$exonType, levels = unique(significantDomainsOR$exonType))
levels(significantDomainsOR$exonType) <- c("EN", "EP")

ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
			axis.title = element_text(size = 10, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 2),
			panel.grid.major = element_blank(), #panel.grid.minor.x = element_line(colour="black", size=1)
			strip.text.x = element_text(size = 10, colour = "brown")
			)

							
ggplot(significantDomainsOR, aes(x = exonType, y = domain)) + 
            geom_point(aes(size = OR, shape = significance)) + scale_shape_manual(values=c(21, 16)) +
			scale_size_continuous(range = c(1,6), breaks = c(1, 5, 10, 20, 50, 100)) +
			geom_hline(aes(yintercept = yLines), color = "grey") + geom_vline(aes(xintercept = xLines), color = "grey") +
    		facet_wrap(~Tissue, scale = "free_x") + theme_bw() + ggTheme


ggsave(file = glue("plotForDomainEnrichmentOfExons.svg"), width =  7, height = 6) 



#orMatrix <- lapply(domainEnrichment, function(x) x[x$domain %in% significantDomainNames, c('domain', 'OR')])
#orMatrix <- do.call("cbind", orMatrix) 
#indexes <- seq(2,ncol(orMatrix), by = 2)
#orMatrix <- orMatrix[ ,indexes] %>% set_rownames((orMatrix[ ,1]))
#orMatrix = sapply(orMatrix, as.numeric) %>% as.data.frame() %>% set_rownames(rownames(orMatrix))
#orMatrix <- do.call(data.frame, lapply(orMatrix, function(x) replace(x, is.infinite(x), 50))) %>% set_rownames((rownames(orMatrix)))
              


#orMatrix <- significantDomainsOR %>% subset(., Tissue == "Liver", select = c(exonType, domain, OR)) %>%
#    reshape2::dcast(data = .,formula = domain ~ exonType, value.var = "OR")
    #tidyr::spread(., key = domain, value = OR)





orMatrix <- significantDomainsOR %>% #subset(., Tissue == "Liver", select = c(exonType, domain, OR)) %>%
    reshape2::dcast(data = .,formula = domain + Tissue ~ exonType, value.var = "OR")

fdrMatrix <- significantDomainsOR %>% #subset(., Tissue == "Liver", select = c(exonType, domain, OR)) %>%
    reshape2::dcast(data = .,formula = domain + Tissue ~ exonType, value.var = "FDR")

orMatrix[is.infinite(orMatrix$embPos), 'embPos'] <- 150 #max(orMatrix$embPos[is.finite(orMatrix$embPos)])
orMatrix[is.infinite(orMatrix$embNeg), 'embNeg'] <- 150 #max(orMatrix$embNeg[is.finite(orMatrix$embNeg)])

significance <- ifelse(fdrMatrix$embNeg < 0.1 | fdrMatrix$embPos < 0.1, "significant", "not significant")
orMatrix$significance <- significance

ggplot(data = orMatrix, aes(x = embPos, y = embNeg)) + geom_point() + facet_wrap(~Tissue, scale = "free") + 
		coord_cartesian(xlim = c(-50,180), ylim = c(0,180)) + #geom_label_repel(aes(label = domain), size = 2) + 
		geom_label_repel(data = subset(orMatrix, orMatrix$significance == "significant"), aes(label = subset(orMatrix, orMatrix$significance == "significant")$domain),
		    	size = 2, fontface = 'bold', max.overlaps = 50, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines" #nudge_y = 0.1,
			))

