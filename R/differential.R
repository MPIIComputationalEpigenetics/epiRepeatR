rankDESeqTab <- function(dm){
	if (is.element("DESeqResults", class(dm))){
		dm <- data.frame(dm)
	}
	rankMat <- cbind(
		rank(-abs(dm[,"log2FoldChange"]), na.last="keep", ties.method="min"),
		rank(dm[,"pvalue"], na.last="keep", ties.method="min")
	)
	dm[,"cRank"] <- rowMaxs(rankMat, na.rm=TRUE)
	dm[dm[,"cRank"]==-Inf,"cRank"] <- NA
	dm[,"cRank_rerank"] <- rank(dm[,"cRank"], na.last="keep", ties.method="min")
	return(dm)
}

computeDiffMat.counts.deseq2 <- function(countTab, sampleAnnot, designF, colName, g1Name, g2Name){
	require(DESeq2)
	dds <- DESeqDataSetFromMatrix(countData=countTab, colData=sampleAnnot, design=designF)
	dds <- DESeq(dds)
	diffRes <- results(dds, contrast=c(colName, g1Name, g2Name))
	dm <- rankDESeqTab(diffRes)
	normCounts <- counts(dds, normalized=TRUE)
	dm[,"baseMean.g1"] <- rowMeans(normCounts[,sampleAnnot[,colName]==g1Name], na.rm=TRUE)
	dm[,"baseMean.g2"] <- rowMeans(normCounts[,sampleAnnot[,colName]==g2Name], na.rm=TRUE)
	dm[,"log2BaseMean"] <- log2(dm[,"baseMean"])
	dm[dm[,"log2BaseMean"]==-Inf,"log2BaseMean"] <- NA
	dm <- dm[,c("log2BaseMean", "log2FoldChange", "pvalue", "padj", "cRank", "cRank_rerank",  "baseMean.g1", "baseMean.g2", "baseMean", "lfcSE", "stat")]
	#head(cbind(diffScores$baseMean, rowMeans(diffScores[,c("baseMean.g1","baseMean.g2")]), diffScores$baseMean.g1, diffScores$baseMean.g2, diffScores$baseMean*(2^diffScores$log2FoldChange)))
	#head(cbind(log2(diffScores$baseMean.g1/diffScores$baseMean.g2), diffScores$log2FoldChange))
	return(dm)
}
