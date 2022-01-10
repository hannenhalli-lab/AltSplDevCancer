library(preprocessCore)
exonGeneNames <- function(exonList, GenesDf) {
	exonGenes <-  exonList %>% as.character() %>% gsub("\\..+", "", .)
	indexes <- match(exonGenes, GenesDf[ ,'GeneID'])
	exonGenes <- GenesDf[indexes, c('GeneName')] %>% as.character()
	geneType <- GenesDf[indexes, c('GeneType')] %>% as.character()
	return(cbind(Exon = exonList, GeneName = exonGenes, geneType = geneType))
}


makeExonBed <- function(InputFile, GenesDf) {
	#Genes <- gsub("\\..+", "", InputFile$Exon)
	Genes <- exonGeneNames(exonList = InputFile$Exon, GenesDf = GenesDf) %>% data.frame()
	strand <- InputFile$Exon %>% gsub(".+:", "", .)
	Exons <- InputFile$Exon %>% 
		gsub(".+SE:|[:-]\\d+:[+-]", "", .) %>%
		gsub("\\d+-", "", .) %>%
		strsplit(., ":") %>% unlist() %>% matrix(., ncol = 3, byrow = T) %>% 
		data.frame() %>% mutate(exon = InputFile$Exon, geneName = Genes$exonGenes, strand = strand)
	names(Exons)[1:3] <- c("chr", "start", "end")
	Exons$chr <- as.character(Exons$chr)
	Exons$start <- Exons$start %>% as.character %>% as.numeric()
	Exons$end <- Exons$end %>% as.character %>% as.numeric()
	#Exons <- Exons[order(Exons$chr, Exons$start), ]
	Exons$start <- Exons$start - 1
	return(Exons)
}


processEpxressionFile <- function(dataFrame, Transcripts, replicates, normalize = T) {
  if (normalize) {
    print("performing quantile normalization")
  } else {
    print("processing without quantile normalization")
  }
	rownames(dataFrame) <- gsub("\\.\\d+", "", rownames(dataFrame))
	matches <- match(rownames(dataFrame), Transcripts$tx_id)
	Genes <- Transcripts[matches, 'gene_name']
	dataFrame <- cbind("GeneName" = Genes, dataFrame)
	GeneLevelTPM <- aggregate(.~GeneName, data = dataFrame, sum, na.rm = T, na.action = na.pass) 
	if (normalize) {
	  GeneLevelTPM <- GeneLevelTPM[ ,-1] %>% as.matrix() %>% normalize.quantiles() %>% 
	    set_rownames(GeneLevelTPM$GeneName) %>% 
	    set_colnames(colnames(GeneLevelTPM)[-1]) %>% data.frame(check.names = F)
	}else {
	  GeneLevelTPM <- GeneLevelTPM[ ,-1] %>% set_rownames(GeneLevelTPM$GeneName)
	}
	if (replicates) {
		GeneLevelTPM <- GeneLevelTPM %>% t() %>% data.frame() %>%
				mutate(Stage = gsub("_Rep\\d+", "", colnames(GeneLevelTPM))) %>% 
				aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
				set_rownames(.$Stage) %>% .[,-1] %>% t() %>% data.frame() %>%
				set_rownames(rownames(GeneLevelTPM))
	}	
	
	#GeneLevelTPM <- GeneLevelTPM[ ,mixedsort(names(GeneLevelTPM))]
	return(GeneLevelTPM)
}

processPSIfile <- function(psiFile) {
	data.frame(t(psiFile)) %>% 
	mutate(Stage = gsub("_Rep\\d+", "", names(psiFile))) %>% 
	aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
	set_rownames(.$Stage) %>% 
	.[,-1] %>% t() %>% data.frame() %>%
	set_rownames(rownames(psiFile))
}

naFilter <- function(dataFrame, cutoff) {
	#nsize <- floor(ncol(dataFrame) * cutoff)
  nsize <- ceiling(ncol(dataFrame) * cutoff)
	indexes <- apply(dataFrame, 1, function(x) sum(is.na(x))) < nsize
	dataFrame[indexes, ]
}

