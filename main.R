library(hgu133plus2.db)
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(readxl)
library(writexl)
library(tibble)
library(scales)
library(viridis)
library(tidyr)

rawCountsFile <- "RawCounts.xlsx"
sampleDataFile <- "SampleData.xlsx"
Compare <- "Disease"
threshold <- 5
saveFolder <- "Results"

dir.create(saveFolder)

k <- keys(hgu133plus2.db,keytype="PROBEID")
mapping<-AnnotationDbi::select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="PROBEID")
mapping<-mapping[!duplicated(mapping$PROBEID), ]
rownames(mapping)<-mapping$PROBEID
mapping <- subset(mapping, select = -PROBEID)

rawCounts <- data.frame(read_excel(rawCountsFile))
sampleData <- data.frame(read_excel(sampleDataFile))

sampleData <- as.data.frame(apply(sampleData, 2, function(x) gsub("\\s+", "_", x)))

ProbeID <- rawCounts$ProbeID
rawCounts <- as.matrix(subset(rawCounts, select = -ProbeID))
rownames(rawCounts) <- ProbeID

SampleID = sampleData$SampleID
rownames(sampleData) <- SampleID

sampleData <- subset(sampleData, select = -SampleID)

sampleData[[Compare]] <- factor(sampleData[[Compare]])

rawCounts <- rawCounts[,unique(rownames(sampleData))]
print(all(colnames(rawCounts) == rownames(sampleData)))

sampleData[[Compare]] <- factor(sampleData[[Compare]], levels=unique(sampleData[[Compare]]))

exp_raw <- log2(rawCounts)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame("PC1" = PCA_raw$x[,1], 
                     "PC2" = PCA_raw$x[,2],
                     "tmp" = sampleData[,Compare])

colnames(dataGG) = c("PC1","PC2",Compare)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes_string(shape = Compare)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15))

ggsave(file.path(saveFolder, "PCA.pdf"))

ggplot(stack(data.frame(exp_raw)), aes(x = ind, y = values, fill=ind)) +
  xlab("") +ylab("")  + theme(legend.position = "none") + geom_boxplot(outlier.size = 0.1) +  
  labs(title="Boxplot of log2-intensitites for the raw data") 

ggsave(file.path(saveFolder, "boxplot_log2_intensities.pdf"))

my_deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ Disease)
my_deseq2Data <- my_deseq2Data[rowSums(counts(my_deseq2Data)) > threshold, ]
my_deseq2Data <- DESeq(my_deseq2Data)

normalized_counts<-log2(counts(my_deseq2Data, normalized=TRUE))

row_medians_assayData <- Biobase::rowMedians(as.matrix(normalized_counts))

RLE_data <- sweep(normalized_counts, 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, sample_array, log_expr_dev)

ggplot2::ggplot(RLE_data_gathered, aes(x = sample_array, y = log_expr_dev, fill=sample_array)) + 
  xlab("") +ylab("")  + theme(legend.position = "none") +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = quantile(RLE_data_gathered$log_expr_dev, c(0.05, 0.95))) + 
  labs(title="Boxplot of log2-expression deviation")

ggsave(file.path(saveFolder, "boxplot_log2_expression_deviation.pdf"))

exp <- counts(my_deseq2Data, normalized=TRUE)
PCA <- prcomp(t(exp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame("PC1" = PCA_raw$x[,1], 
                     "PC2" = PCA_raw$x[,2],
                     "tmp" = sampleData[,Compare])

colnames(dataGG) = c("PC1","PC2",Compare)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes_string(shape = Compare)) +
  ggtitle("PCA plot of the normalized expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15))

ggsave(file.path(saveFolder, "PCA_normalized.pdf"))

annotation_for_heatmap <- data.frame("tmp" = sampleData[,Compare])
colnames(annotation_for_heatmap) = c(Compare)
row.names(annotation_for_heatmap) <- row.names(sampleData)

dists <- as.matrix(dist(t(exp), method = "manhattan"))

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
         main = "Clustering heatmap",
         filename=file.path(saveFolder, "heatmap.pdf"))

medians <- rowMedians(log2(exp))

pdf(file.path(saveFolder, "hist_median_intensities.pdf"))
hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = threshold, col = "coral4", lwd = 3)
dev.off()

my_deseq2Results <- results(my_deseq2Data, contrast=c(Compare, levels(sampleData[[Compare]])))

deseq2ResDF <- as.data.frame(my_deseq2Results)
deseq2ResDF <- add_column(deseq2ResDF, SYMBOL = mapping[rownames(deseq2ResDF),], .after = 0)
deseq2ResDF <- add_column(deseq2ResDF, PROBEID = rownames(deseq2ResDF), .after = 0)
write_xlsx(deseq2ResDF, file.path(saveFolder, "DiffExprRes.xlsx"))

deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < 0.1, "Significant", NA)

ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) +
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw()

ggsave(file.path(saveFolder, "fold_vs_count.pdf"))

ggplot(deseq2ResDF, aes(log2FoldChange, baseMean, colour=significant)) +
  geom_point(size=1) + scale_x_continuous(limits=c(-3, 3), oob=squish) +
  scale_y_log10() + geom_vline(xintercept = 0, colour="tomato1", size=2) +
  labs(y="mean of normalized counts", x="log2 fold change") +
  scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

ggplot(deseq2ResDF, aes(log2FoldChange, baseMean, colour=padj)) +
  geom_point(size=1) + scale_x_continuous(limits=c(-3, 3), oob=squish) +
  scale_y_log10() + geom_vline(xintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(y="mean of normalized counts", x="log2 fold change") +
  scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw()

ggsave(file.path(saveFolder, "count_vs_fold.pdf"))

logdeseq2ResDF = deseq2ResDF
logdeseq2ResDF$padj = -log10(logdeseq2ResDF$padj)

ggplot(logdeseq2ResDF, aes(log2FoldChange, padj)) +
  geom_point(size=1) + scale_x_continuous(limits=c(-3, 3), oob=squish) +
  labs(y="-log10 P-value", x="log2 fold change")

ggsave(file.path(saveFolder, "volcano.pdf"))