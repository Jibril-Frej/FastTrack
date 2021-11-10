# Creates and returns a mapping between RNA_ID and SYMBOL 
build_mapping <- function(ID_TYPE){
  k <- keys(hgu133plus2.db,keytype=ID_TYPE)
  mapping<-AnnotationDbi::select(hgu133plus2.db, keys=k, columns=c("GENENAME"), keytype=ID_TYPE)
  mapping<-mapping[!duplicated(mapping[[ID_TYPE]]), ]
  rownames(mapping) <- mapping[[ID_TYPE]]
  mapping_G <- subset(mapping, select = -get(ID_TYPE))
  mapping<-AnnotationDbi::select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype=ID_TYPE)
  mapping<-mapping[!duplicated(mapping[[ID_TYPE]]), ]
  rownames(mapping) <- mapping[[ID_TYPE]]
  mapping_S <- subset(mapping, select = -get(ID_TYPE))
  return(list(mapping_G = mapping_G, mapping_S = mapping_S))
}

# Read the xlsx files containing the raw counts and sample information and returns two objects taht DESEQ2 can process 
read_raw_data <- function(rawCountsFile, sampleDataFile){
  
  # Read the xlsx files into data.frames
  rawCounts <- data.frame(read_excel(rawCountsFile))
  sampleData <- data.frame(read_excel(sampleDataFile))
  
  # Changes rawCounts into a matrix for DESEQ2 and use RNA_ID as rownames
  RNA_ID <- rawCounts$RNA_ID
  rawCounts <- as.matrix(subset(rawCounts, select = -RNA_ID))
  rownames(rawCounts) <- RNA_ID
  
  
  # Replaces all spaces in sampleData with underscore (needed for DESEQ2)
  sampleData <- as.data.frame(apply(sampleData, 2, function(x) gsub("\\s+", "_", x)))
  
  # Use SampleID as rownames
  SampleID = sampleData$SampleID
  rownames(sampleData) <- SampleID
  sampleData <- subset(sampleData, select = -SampleID)
  
  sampleData[[Compare]] <- factor(sampleData[[Compare]])
  
  # Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
  rawCounts <- rawCounts[,unique(rownames(sampleData))]
  print(all(colnames(rawCounts) == rownames(sampleData)))
  
  sampleData[[Compare]] <- factor(sampleData[[Compare]], levels=unique(sampleData[[Compare]]))
  
  return(list(rawCounts = rawCounts, sampleData = sampleData))
  
}


# Generates and saves a pdf of the PCA plot of the log-transformed raw expression data
PCA <- function(logRawCounts, sampleData, Compare, saveFolder){
  
  PCA_raw <- prcomp(t(logRawCounts), scale. = FALSE)
  
  percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  dataGG <- data.frame("PC1" = PCA_raw$x[,1], 
                       "PC2" = PCA_raw$x[,2],
                       "tmp" = sampleData[,Compare])
  
  colnames(dataGG) = c("PC1","PC2",Compare)
  
  p <- ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes_string(shape = Compare)) +
    ggtitle("PCA plot of the log2 raw counts") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio) +
    scale_shape_manual(values = c(4,15))
  
  print(p)
  
  ggsave(file.path(saveFolder, "PCA.pdf"))
  
  
}


box_intensities <- function(logRawCounts, saveFolder){
  p <- ggplot(stack(data.frame(logRawCounts)), aes(x = ind, y = values, fill=ind)) +
    xlab("") +ylab("log2 raw counts")  + theme(legend.position = "none") + geom_boxplot(outlier.size = 0.1) +  
    labs(title="Boxplot of the log2 raw counts") +
    theme(axis.text.x = element_text(angle = 90))
  print(p)
  ggsave(file.path(saveFolder, "boxplot_log2_raw_counts.pdf"))
}


make_DEseq2DataSet <-function(rawCounts, sampleData, Compare, threshold, saveFolder){
  
  deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, 
                                          design=as.formula(paste("~", Compare, sep="")))
  deseq2Data <- DESeq(deseq2Data)
  
  normalized_counts <- log2(counts(deseq2Data, normalized=TRUE) +1)
  
  medians <- rowMedians(normalized_counts)
  
  #pdf(file.path(saveFolder, "hist_median_intensities.pdf"))
  hist_res <- hist(medians, 100, col = "cornsilk", freq = TRUE,
                   main = "Histogram of the log2 normalized counts",
                   border = "antiquewhite4",
                   xlab = "Median normalized raw counts",
                   ylab = "Number of RNA strands")
  abline(v = log2(threshold), col = "coral4", lwd = 3)
  #plot(hist_res)
  #dev.off()
  dev.print(pdf, file.path(saveFolder, 'hist_median_intensities.pdf'))

  deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, 
                                       design=as.formula(paste("~", Compare, sep="")))
  
  deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > threshold, ]
  deseq2Data <- DESeq(deseq2Data)
  return(deseq2Data)
}

box_deviation <- function(deseq2Data, saveFolder){
  
  normalized_counts<-log2(counts(deseq2Data, normalized=TRUE)+1)
  
  row_medians_assayData <- Biobase::rowMedians(as.matrix(normalized_counts))
  
  RLE_data <- sweep(normalized_counts, 1, row_medians_assayData)
  RLE_data <- as.data.frame(RLE_data)
  RLE_data_gathered <- tidyr::gather(RLE_data, sample_array, log_expr_dev)
  
  p <- ggplot(RLE_data_gathered, aes(x = sample_array, y = log_expr_dev, fill=sample_array)) + 
    xlab("") +ylab("deviation")  + theme(legend.position = "none") +
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(RLE_data_gathered$log_expr_dev, c(0.05, 0.95))) + 
    labs(title="Boxplot of log2 normalized counts deviation")+
    theme(axis.text.x = element_text(angle = 90))
  
  print(p)
  
  ggsave(file.path(saveFolder, "boxplot_normalized_deviation.pdf"))
  
}

normalized_PCA <- function(normalized_counts, sampleData, Compare, saveFolder){
  
  PCA <- prcomp(t(normalized_counts), scale = FALSE)
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  dataGG <- data.frame("PC1" = PCA$x[,1], 
                       "PC2" = PCA$x[,2],
                       "tmp" = sampleData[,Compare])
  
  colnames(dataGG) = c("PC1","PC2",Compare)
  
  p <- ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes_string(shape = Compare)) +
    ggtitle("PCA plot of the log2 normalized counts") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = sd_ratio) +
    scale_shape_manual(values = c(4,15))
  
  print(p)
  
  ggsave(file.path(saveFolder, "PCA_normalized.pdf"))
  
}

heatmap <- function(normalized_counts, sampleData, Compare, saveFolder){
  
  annotation_for_heatmap <- data.frame("tmp" = sampleData[,Compare])
  colnames(annotation_for_heatmap) = c(Compare)
  row.names(annotation_for_heatmap) <- row.names(sampleData)
  
  dists <- as.matrix(dist(t(normalized_counts), method = "manhattan"))
  
  rownames(dists) <- row.names(sampleData)
  hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
  colnames(dists) <- NULL
  diag(dists) <- NA
  
  pheatmap(dists, col = (hmcol), 
           annotation_row = annotation_for_heatmap,
           legend = TRUE, 
           treeheight_row = 0,
           legend_breaks = c(min(dists, na.rm = TRUE), 
                             max(dists, na.rm = TRUE)),
           legend_labels = (c("small distance", "large distance")),
           main = "Clustering heatmap")
  
  dev.print(pdf, file.path(saveFolder, 'heatmap.pdf'))

}


diff_Analysis <- function(deseq2Results, deseq2Data, sampleData, Compare, mapping, saveFolder){
  
  deseq2Results <- results(deseq2Data, contrast=c(Compare, levels(sampleData[[Compare]])))
  deseq2ResDF <- as.data.frame(deseq2Results)
  deseq2ResDF <- add_column(deseq2ResDF, SYMBOL = mapping$mapping_S[rownames(deseq2ResDF),], .after = 0)
  deseq2ResDF <- add_column(deseq2ResDF, GENENAME = mapping$mapping_G[rownames(deseq2ResDF),], .after = 0)
  deseq2ResDF <- add_column(deseq2ResDF, ID = rownames(deseq2ResDF), .after = 0)
  write_xlsx(deseq2ResDF, file.path(saveFolder, "DiffExprRes.xlsx"))
  deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < 0.1, "Significant", NA)
  
  return(deseq2ResDF)
  
}


fold_vs_count <- function(deseq2ResDF, saveFolder){
  
  ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + 
    geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) +
    scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) +
    labs(x="mean of normalized counts", y="log fold change") +
    scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()
  
  p<-ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) +
    geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) +
    scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
    labs(x="mean of normalized counts", y="log fold change") +
    scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw()
  print(p)
  ggsave(file.path(saveFolder, "fold_vs_count.pdf"))
  
}

count_vs_fold <- function(deseq2ResDF, saveFolder){
  
  ggplot(deseq2ResDF, aes(log2FoldChange, baseMean, colour=significant)) +
    geom_point(size=1) + scale_x_continuous(limits=c(-3, 3), oob=squish) +
    scale_y_log10() + geom_vline(xintercept = 0, colour="tomato1", size=2) +
    labs(y="mean of normalized counts", x="log2 fold change") +
    scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()
  
  p<-ggplot(deseq2ResDF, aes(log2FoldChange, baseMean, colour=padj)) +
    geom_point(size=1) + scale_x_continuous(limits=c(-3, 3), oob=squish) +
    scale_y_log10() + geom_vline(xintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
    labs(y="mean of normalized counts", x="log2 fold change") +
    scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw()
  
  print(p)
  ggsave(file.path(saveFolder, "count_vs_fold.pdf"))


}

volcano <- function(deseq2ResDF, saveFolder){

  logdeseq2ResDF = deseq2ResDF
  logdeseq2ResDF$padj = -log10(logdeseq2ResDF$padj)
  
  p<-ggplot(logdeseq2ResDF, aes(log2FoldChange, padj)) +
    geom_point(size=1) + scale_x_continuous(limits=c(-3, 3), oob=squish) +
    labs(y="-log10 P-value", x="log2 fold change")
  
  print(p)
  
  ggsave(file.path(saveFolder, "volcano.pdf"))


}