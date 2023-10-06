#!/usr/bin/env -S Rscript

# Author: Yutathkarn Coles
# Description: Plot gene counts of specific genes

# Variables
META <- "meta-ptr.tsv"
MATRIX <- "ptr.tsv"
GENES <- c("data/mco-genes.txt", "data/effector-genes.txt")
CONTRAST <- "groupCultivar"

OUT <- "out"
THREADS <- 4

# Install packages using autoload.R
source("bin/autoload.R")

# Set parameters for parallelisation, requires doParallel and BiocParallel
if (!exists("THREADS")) {
  THREADS <- detectCores()
}
registerDoParallel(cores = THREADS)
register(MulticoreParam(workers = THREADS))

# Load command variables
argv <- commandArgs(trailingOnly = TRUE)

# Load settings
## conf <- read.config(file = "settings.toml")

# Load data
source("bin/common.R")
source("bin/metadataFungi.R") # FIXME: Lost file to corruption
source("bin/functions.R")

matrix <- read.csv(MATRIX, sep = "\t", header = TRUE)
meta <- read.csv(META, sep = "\t", header = TRUE)
normGene <- "actin"

samplesNumber <- which(meta[, "control"] == "disease")
# TODO: Find more portable format
groups <- meta[, CONTRAST][samplesNumber][replicates * (1:(length(samplesNumber) / replicates))]
groupsOrdered <- as.vector(sapply(
  cultivarsOrdered,
  function(x) {
    return(c(paste0(x, "Early"), paste0(x, "Late")))
  }
))
samples <- c()
for (i in 1:length(groups)) {
  samples[3 * i - 2] <- paste0(groups[i], samplesNumber[3 * i - 2])
  samples[3 * i - 1] <- paste0(groups[i], samplesNumber[3 * i - 1])
  samples[3 * i - 0] <- paste0(groups[i], samplesNumber[3 * i - 0])
}
samplesOrdered <- c()
for (i in 1:length(groupsOrdered)) {
  l <- nchar(groupsOrdered[i])
  for (j in 1:length(samples)) {
    x <- substr(samples[j], 1, l)
    if (x == groupsOrdered[i]) {
      samplesOrdered <- append(samplesOrdered, samples[j])
    }
  }
}

genes <- c()
ncat <- list()
buffer <- c()
for (i in 1:length(GENES)) {
  cat <- strsplit(basename(GENES[i]), split = "[.]")[[1]][1]
  buffer <- as.vector(fread(GENES[i], header = FALSE))
  genes <- append(genes, buffer[[1]])
  print(genes)
  ncat[[i]] <- list()
  ncat[[i]][[1]] <- length(buffer[[1]])
  ncat[[i]][[2]] <- cat
}
n <- 1
for (i in 1:length(ncat)) {
  print(n)
  names(genes)[n:(n + ncat[[i]][[1]] - 1)] <- ncat[[i]][[2]]
  n <- n + ncat[[i]][[1]]
}

genes <- appendGene("PtrM4_092750", "actin")
genes <- appendGene("PtrM4_118660", "toxa")
normGene <- "actin"

# Remove controls from samples
matrixOrig <- matrix
metaOrig <- meta
matrix <- matrix[, samplesNumber]
meta <- meta[samplesNumber, ]
colnames(matrix) <- rownames(meta)

dds <- loadExperimentData(matrix, meta)
counts <- counts(dds, normalize = TRUE)
rawCounts <- counts(dds, normalize = FALSE)

rlog(dds, blind = FALSE)

# Plotting genes
totalCounts <- allCounts(rawCounts)
totalCountsNorm <- allCounts(counts)

heatmapAllData <- plotGenesHeatmap(counts = rawCounts, totalCounts = totalCounts, mode = "sample")
heatmapAvgData <- plotGenesHeatmap(counts = rawCounts, totalCounts = totalCounts, mode = "group")
heatmapAllDataNorm <- plotGenesHeatmap(counts = counts, totalCounts = totalCountsNorm, mode = "sample")
heatmapAvgDataNorm <- plotGenesHeatmap(counts = counts, totalCounts = totalCountsNorm, mode = "group")

barplotGeneCounts <- heatmapAvgData[[4]]
barplotGeneCountsNorm <- heatmapAvgDataNorm[[4]]

boxplotActin <- plotGeneBox(gene = genes["actin"], counts = counts)
boxplotToxA <- plotGeneBox(gene = genes["toxa"], counts = counts)
boxplotNormToxA <- plotGeneNormBox(gene = genes["toxa"], counts = counts, norm = genes[normGene])

grid.arrange(heatmapAvgData[[1]])
grid.arrange(heatmapAvgDataNorm[[1]])

g <- grid.arrange(
  rbind(
    ggplotGrob(barplotGeneCounts),
    ggplotGrob(boxplotActin),
    ggplotGrob(boxplotToxA),
    ggplotGrob(boxplotNormToxA)
  ),
  top = "Raw counts, normalised expression and ratios"
)

# PCA plots without edgeR or DESeq2
pca0 <- PCA(rawCounts, graph = FALSE)
genesForPCA <- names(
  head(sort(abs(pca0[["ind"]][[1]][, 1]), decreasing = TRUE), n = 2000)
)

## ## pca1 <- PCA(matrix[genes, ], graph = FALSE)
## pca1 <- PCA(matrix[genesForPCA, ], graph = FALSE)
## pca1.pc1 <- pca1[["ind"]][["coord"]][, 1]
## pca1.pc2 <- pca1[["ind"]][["coord"]][, 2]
## pca1.var <- (pca1[["var"]][["coord"]] %>% data.frame())[, 1:2]
## colnames(pca1.var) <- c("pc1", "pc2")
## pca1.var[, "sample"] <- samples
## pca1.var[, "cultivar"] <- rep(cultivars, each = 2 * replicates)
## pca1.var[, "group"] <- rep(groups, each = replicates)
## pca1.var[, "time"] <- meta[, "time"]

## p <- ggplot(data = pca1.var, aes(x = pc1, y = pc2, color = cultivar, shape = time)) +
##   geom_point(size = 3)
## p

## g1 <- ggplotGrob(heatmapAvgData[[4]])
## g2 <- ggplotGrob(barplotActin)
## g <- rbind(g1, g2)
## grid.arrange(g)

# Average counts per group
## avgCounts <- averageCounts(counts)
## plotGenesAllHeatmap(allCounts = allCounts, counts = counts)

# Differential expression analysis
## contrasts <- getContrastsFungi(dds)
## sdeg <- getSDEG(contrasts)

## plots <- plotContrastsHeatmap(contrasts)
## saveContrastsHeatmap(plots)

## # Do pairwise comparison of significant genes from each contrast
## x <- sdeg
## y <- list()
## for (i in 1:length(x)) {
##   for (j in i:length(x)) {
##     g <- intersect(rownames(x[[i]]), rownames(x[[j]]))
##   }
## }
