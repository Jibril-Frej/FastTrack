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

source("utils.R")

# Path of the xlxs file containing the raw non-normalized sequence read counts at either the gene or transcript level
rawCountsFile <- "RawCounts.xlsx"

# Path of the xlxs file containing  the properties of the samples
sampleDataFile <- "SampleData.xlsx"

# Indicates which sample property to compare
Compare <- "Disease"

# Threshold used to discard the genes/transcript that have a small amount of read count
threshold <- 30

# Path of the folder where the results and plots will be saved
saveFolder <- "Results"
dir.create(saveFolder)

mapping <- build_mapping()

rawData <- read_raw_data(rawCountsFile, sampleDataFile)
logRawCounts <- log2(rawData$rawCounts)


PCA(logRawCounts, rawData$sampleData, Compare, saveFolder)

box_intensities(logRawCounts, saveFolder)

deseq2Data <- make_DEseq2DataSet(rawData$rawCounts, rawData$sampleData, Compare, threshold)

box_deviation(deseq2Data, saveFolder)

normalized_counts <- log2(counts(deseq2Data, normalized=TRUE))

normalized_PCA(normalized_counts, rawData$sampleData, Compare, saveFolder)

heatmap(normalized_counts, rawData$sampleData, Compare, saveFolder)

deseq2ResDF <- diff_Analysis(deseq2Results, deseq2Data, rawData$sampleData, Compare, mapping, saveFolder)


fold_vs_count(deseq2ResDF, saveFolder)
count_vs_fold(deseq2ResDF, saveFolder)
volcano(deseq2ResDF, saveFolder)