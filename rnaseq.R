# Author: Yutathkarn Coles
# Date: Wed Aug 16 14:26:05 2023

# Wheat and fungi RNA-seq analysis

# File handling with base R
source("lib/path.R")

# Install packages using autoload.R
source(formatPath("lib", "autoload.R"))

# Load settings
source(formatPath("lib", "settings.R"))

# Load globals
for (i in list.files(formatPath("lib", "globals"))) {
  source(formatPath("lib", "globals", i))
}

source(formatPath("lib", "parallel.R"))
source(formatPath("lib", "classes.R"))
source(formatPath("lib", "colours.R"))

source(formatPath("bin", "data.R"))
source(formatPath("bin", "metadata.R"))
source(formatPath("bin", "deseq.R"))

for (i in list.files(formatPath("bin", "graphics"))) {
  source(formatPath("bin", "graphics", i))
}

# Very important for bar plots!
set.seed(1)

geneList <- loadGeneList(genes)
addSelectGenes(geneList, genesSelect)
expData <- loadData(mat[, samples], meta[samples, ], contrast, factors, nreps, species)
expFactors <- getFactors(factorList, factors)
setFactorSamples(expData, expFactors, contrast)
setFactor(expData, expFactors, list(hkg = hkg, nreps = expData[["nreps"]], factors = factors))

fixData(expData)
dds <- loadExperimentData(expData, geneList, threshold = countThreshold)

# Counts
countsRaw <- counts(dds, normalized = FALSE)
countsNorm <- counts(dds, normalized = TRUE)
countsRL <- assay(rlog(dds))
countsVS <- assay(vst(dds))
colnames(countsRL) <- meta[meta[, "control"] == "disease", "sample"]
colnames(countsVS) <- meta[meta[, "control"] == "disease", "sample"]

# Counts per column (sample)
totalCountsAll <- getTotalCounts(countsRaw, 1)
totalCountsAllNorm <- getTotalCounts(countsNorm, 1)
totalCountsAvg <- getTotalCounts(countsRaw, expFactors[["factorList"]][["nreps"]])
totalCountsAvgNorm <- getTotalCounts(countsNorm, expFactors[["factorList"]][["nreps"]])

# Plots
plotDefs <- getPlotDefs(confPlots[["plots"]])
pTotalCounts <- plotTotalCounts(totalCountsAll, expFactors, plotDefs)
pTotalCountsBars <- plotTotalCountsBars(totalCountsAll, expFactors, plotDefs)
pTotalCountsNorm <- plotTotalCounts(totalCountsAllNorm, expFactors, plotDefs)
pTotalCountsBarsNorm <- plotTotalCountsBars(totalCountsAllNorm, expFactors, plotDefs)
# TODO: Change heatmaps to have boxes around gene groups
pCountsHeatmap <- plotCountsHeatmap(countsRaw, geneList, expFactors, plotDefs)

for (i in 1:length(geneList[["genes"]])) {
  print(i)
  countsRaw[geneList[["genes"]][i], ]
}

## pCountsHeatmapAll <- plotCountsHeatmap(countsRaw, geneList, expFactors, plotDefs, subset = FALSE, threshold = 1000)
pCountsHeatmapAll50 <- plotCountsHeatmapAll(countsRaw, expFactors, plotDefs, threshold = 50)

p1 <- grid.arrange(pTotalCounts, pTotalCountsBars, nrow = 1)
p2 <- grid.arrange(pTotalCountsNorm, pTotalCountsBarsNorm, nrow = 1)

savePlotGrid(p1, formatPath(expData[["species"]], "pTotalCounts_BOTH.png", root = FALSE), plotDefs)
savePlotGrid(p2, formatPath(expData[["species"]], "pTotalCountsNorm_BOTH.png", root = FALSE), plotDefs)
saveHeatmap(pCountsHeatmapAll50, formatPath(expData[["species"]], "pCountsHeatmapAll_top50.png", root = FALSE), plotDefs, width = 1000, height = 1000)
savePlotStd(pCountsHeatmap, formatPath(expData[["species"]], "pCountsHeatmap.png", root = FALSE), plotDefs, width = 3000, height = 2000)

# Specific genes
# Correlation: library size against actin expression
hkgCounts <- countsRaw[geneList[["genes"]][names(geneList[["genes"]]) == hkg], ]
pPairCountsScatter <- plotPairCountsScatter(hkgCounts, totalCountsAll, expFactors, plotDefs)
hkgCountsNorm <- countsRaw[geneList[["genes"]][names(geneList[["genes"]]) == hkg], ]
pPairCountsScatterNorm <- plotPairCountsScatter(hkgCountsNorm, totalCountsAllNorm, expFactors, plotDefs)

p3 <- grid.arrange(pPairCountsScatter, pPairCountsScatterNorm, nrow = 1)
savePlotGrid(p3, formatPath(expData[["species"]], "pPairCountsScatter_BOTH.png", root = FALSE), plotDefs, width = 3000, height = 2000)


# Venn diagrams
results <- getResultsList(dds, expData, expFactors)
sdegList <- getDifferentialGenes(results, expFactors)
plotListAll <- plotDifferentialGenes(sdegList, expFactors, plotDefs, offsetR = 0, offsetC = 0, debug = FALSE)
plotList <- plotListAll[[1]]
legendList <- plotListAll[[2]]
saveVennDiagrams(species, plotList, legendList, expFactors, plotDefs)


# Get functional distribution of significant genes


# Expression of top 100 differentially expressed genes from each pair


# Export replicates for Clust
exportToClust(
  formatPath(conf[["data"]][["data"]], paste0(expData[["species"]], "Replicates.txt"), root = FALSE),
  basename(formatPath(fileMat, root = FALSE)),
  contrast, meta[meta[, "control"] == "disease", ]
)
exportToClust(
  formatPath(conf[["data"]][["data"]], paste0(expData[["species"]], "ReplicatesEarly.txt"), root = FALSE),
  basename(formatPath(fileMat, root = FALSE)),
  contrast, meta[meta[, "control"] == "disease" & meta[, "time"] == "early", ]
)
exportToClust(
  formatPath(conf[["data"]][["data"]], paste0(expData[["species"]], "ReplicatesLate.txt"), root = FALSE),
  basename(formatPath(fileMat, root = FALSE)),
  contrast, meta[meta[, "control"] == "disease" & meta[, "time"] == "late", ]
)

write.table(countsRL, formatPath(conf[["data"]][["data"]], paste0(species, "RL.tsv")), sep = "\t", quote = FALSE)
write.table(countsVS, formatPath(conf[["data"]][["data"]], paste0(species, "VS.tsv")), sep = "\t", quote = FALSE)

# TODO: Clusters of log2 changes

# Expression of clusters: clust
# NOTE: Run Clust before running this
# Load clusters (see globals.R for clustRuns)
modules <- getModules(clustRuns)
goMap <- readMappings(fileGoMap)

# Map cluster genes to their annotations using topGO
# Reversed order of iterators... need to reinitialise lists
clustPlots <- list()
clustPlotsTG <- list()

clustPlotsTG <- plotExprProfilesGo(clustPlotsTG, goMap, globalColours = TRUE)
## clustPlots <- plotExprProfilesGoRestricted(clustPlots, goMap, globalColours = FALSE, "aspergillus")

# Aspergillus
# Line graphs for each cluster
# 1. Get subset of genes for each cluster
# 2. Calculate std. dev of expression of genes
# 3. Read in eigengene values for each cluster and cultivar
# 4. Calculate z-score
clustData <- read.table(
  formatPath(conf[["data"]][["data"]], conf[["data"]][["clust"]], species, "late", "Processed_Data", paste0(species, "RL.tsv_processed.tsv")),
  header = TRUE
)
rownames(clustData) <- clustData[, 1]
colnames(clustData)
clustData <- clustData[, -1]

groups <- expFactors[["factorList"]][["group"]]
groups <- groups[2 * 1:(length(groups) / 2)]
groupsOrdered <- expFactors[["factorList"]][["groupOrdered"]]
groupsOrdered <- groupsOrdered[2 * 1:(length(groupsOrdered) / 2)]

plotListME <- list()
for (i in 1:length(modules)) {
  counts <- countsRL
  plotDefsObj <- plotDefs

  genes <- modules[[i]][which(modules[[i]] %in% rownames(counts))]
  genesClust <- clustData[genes, ]

  d <- reshape2::melt(genesClust)
  d[, "gene"] <- rep(genes, length(groups))
  colnames(d) <- c("group", "zscore", "gene")

  d[, "group"] <- factor(d[, "group"], levels = groupsOrdered)

  p <- ggplot(d, aes(x = group, y = zscore, group = gene)) +
    geom_line(size = 0.5, alpha = 0.8) +
    theme_std(plotDefsObj) +
    theme(
      axis.title.x = element_blank(),
    ) +
    # xlab(plotDefs[["defs"]][["names.xlab.cultivar"]]) +
    # ylab(plotDefs[["defs"]][["names.ylab.zscore"]]) +
    ggtitle(paste("Module", i - 1))

  plotListME[[i]] <- p
}

plotListComplex1 <- list()
for (i in 1:length(modules)) {
  plotListComplex1[[i]] <- list()
  for (j in 1:length(ont)) {
    p <- arrangeGrob(clustPlotsTG[[i]][[j]][[2]], clustPlotsTG[[i]][[j]][[1]], plotListME[[i]], ncol = 1, heights = c(0.75, 0.25, 5.0))
    savePlotGrid(
      p,
      formatPath(species, "clust", tolower(ont[j]), paste0("clust_goDist_module", i - 1, ".png"), root = FALSE),
      plotDefs,
      width = 6000,
      height = 4000,
    )
    plotListComplex1[[i]][[j]] <- p
  }
}

plotListComplex2 <- list()
for (i in 1:length(modules)) {
  plotListComplex2[[i]] <- list()
  for (j in 1:(length(ont))) {
    p <- arrangeGrob(clustPlots[[i]][[j]][[2]], clustPlots[[i]][[j]][[1]], plotListME[[i]], ncol = 1, heights = c(0.75, 0.25, 5.0))
    savePlotGrid(
      p,
      formatPath(species, "clust", tolower(ont[j]), paste0("clust_goDist_aspergillus_module", i - 1, ".png"), root = FALSE),
      plotDefs,
      width = 6000,
      height = 4000,
    )
    plotListComplex2[[i]][[j]] <- p
  }
}

plots <- list()
for (i in 1:(length(ont))) {
  for (j in 1:length(modules)) {
    plots[[j]] <- plotListComplex1[[j]][[i]]
  }
  ## p <- do.call(
  ##   arrangeGrob,
  ##   c(
  ##     plots,
  ##     nrow = 2,
  ##     ## top = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.lines.clustExpr.title\"]]} [@{ont[i]}]")),
  ##     ## top = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.lines.clustExpr.title\"]]} [@{ont[i]}]"))
  ##   )
  ## )
  p <- do.call(
    arrangeGrob,
    c(
      plots,
      nrow = 2
      ## top = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.lines.clustExpr.title\"]]}, Aspergillus DB [@{ont[i]}]")),
      ## bottom = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.xlab.cultivar\"]]}"))
    )
  )
  savePlotGrid(
    p,
    formatPath(species, "clust", tolower(ont[i]), paste0("clust_goDist_all", ".png"), root = FALSE),
    plotDefs
  )
}

for (i in 1:(length(ont))) {
  for (j in 1:length(modules)) {
    plots[[j]] <- plotListComplex2[[j]][[i]]
  }
  p <- do.call(
    arrangeGrob,
    c(
      plots,
      nrow = 2
      ## top = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.lines.clustExpr.title\"]]}, Aspergillus DB [@{ont[i]}]")),
      ## bottom = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.xlab.cultivar\"]]}"))
    )
  )
  savePlotGrid(
    p,
    formatPath(species, "clust", tolower(ont[i]), paste0("clust_goDist_aspergillus_all", ".png"), root = FALSE),
    plotDefs
  )
}

## write.table(
##   goNames[which(goNames[, 1] %in% clust1small), ],
##   formatPath(conf[["data"]][["data"]], paste0(species, "-clust1small.tsv")),
##   quote = FALSE,
##   sep = "\t",
##   col.names = FALSE,
##   row.names = FALSE,
## )

## annotCounterNames <- sapply(names(annotCounter), function(x) {paste(strwrap(x, 20), collapse = "\n")})

## d <- data.frame(count = annotCounter, annot = annotCounterNames)
## ggplot(d, aes(x = count, y = annot)) +
##   geom_col()

# WGCNA
## exportTraitsFungi(formatPath(conf[["data"]][["data"]], paste0(expData[["species"]], "Traits.tsv"), root = FALSE), meta[meta[, "control"] == "disease", ], c(factors, extraFactors))
exportTraitsFungi(formatPath(conf[["data"]][["data"]], paste0(expData[["species"]], "Traits.tsv"), root = FALSE), meta[meta[, "control"] == "disease" & meta[, "time"] == "late", ], c(factors, extraFactors))

## exprData <- t(countsRL)
exprData <- t(countsRL[, rownames(meta[meta[, "control"] == "disease" & meta[, "time"] == "late", ])])
## exprData <- fixDataStructure(exprData)
traits <- read.csv(fileTraits, header = TRUE, sep = "\t")
if (colnames(traits)[length(traits)] == "X") {
  traits <- traits[, -c(length(traits))]
}
traits <- traits[, -1]
rownames(traits) <- rownames(exprData)

t <- traits[, c("rating.SVS", "rating.MS", "rating.MR", "rating.R")]

# TODO: Implement proper checkpointing
## x <- list.files(formatPath(conf[["data"]][["data"]], conf[["data"]][["tom"]]))
## if (length(x) > 0) {
##   for (i in x) load(formatPath(conf[["data"]][["data"]], conf[["data"]][["tom"]], i))
## } else {
cor <- WGCNA::cor

# See below for manual heirachial clustering, blockwiseModules (and similar) recommended for large datasets
exprData <- fixDataStructure(exprData)
net <- blockwiseConsensusModules(
  exprData,
  maxBlockSize = 20000,
  power = 6,
  minModulesSize = 30,
  deepSplit = 2,
  pamRespectsDendro = TRUE,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = FALSE,
  verbose = 5,
)

## col = labels2colors(net[["colors"]])
## plotDendroAndColors(
##   net[["dendrograms"]][[1]],
##   col[net[["blockGenes"]][[1]]],
##   "Module colours",
##   groupLabels = "Gene networks",
##   dendroLabels = FALSE,
##   hang = 0.03,
##   addGuide = TRUE,
##   guideHang = 0.05,
## )

## sizeGrWindow(9, 5)

## TODO: Take unaligned reads from wheat and ptr, for Leon's metagenomics project

## Make plots of reads that aligned and did not align
