#!/usr/bin/env -S Rscript

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

# Iterator bounds
lower1 <- 1
upper1 <- length(expFactors[["factorList"]][["group"]]) - 1
lower2 <- lower1 + 1
upper2 <- length(expFactors[["factorList"]][["group"]])

# Counts
countsRaw <- counts(dds, normalized = FALSE)
countsNorm <- counts(dds, normalized = TRUE)
countsRL <- assay(rlog(dds, blind = FALSE))
countsVS <- assay(vst(dds, blind = FALSE))
colnames(countsRaw) <- as.character(expFactors[["factorList"]][["sample"]])
colnames(countsNorm) <- as.character(expFactors[["factorList"]][["sample"]])
colnames(countsRL) <- as.character(expFactors[["factorList"]][["sample"]])
colnames(countsVS) <- as.character(expFactors[["factorList"]][["sample"]])
getCountsAvg()

# Differentially expressed genes
sdeg <- R6Class(NULL, list(
                        names = c(), dn = c(), up = c(), padj = c(),
                        dedup = function() {
                          self$names <- unique(self$names)
                          self$dn <- unique(self$dn)
                          self$up <- unique(self$up)
                        }))$new()
sdegGroups <- getSDEG(sdeg)
sdegTop <- getSDEGTop(sdegGroups, (function() {topN <<- 100})())

# Get number of differentially expressed genes with annotations
lim <- 0
for (i in lower1:upper1) {
  for (j in (lower2 + lim):upper2) {
    ## printf("i: %-3i\tj: %-3i\n", i, j)
    rownames(sdegGroups[[i]][[j - 1 - lim]]) %in% goMap[, 1]
  }
  lim <- lim + 1
}

countsVS["PtrM4_110920", ]
countsVSAvg["PtrM4_110920", ]

rownames(head(sdegTop %>% arrange(padj), topN)) %>%
  as.list %>%
  fwrite(file = formatPath(conf$data$data, qq("sdegTop@{topN}.tsv")), sep = "\t", quote = FALSE)


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
# TODO: Change heatmap to use pheatmap
pCountsHeatmap <- plotCountsHeatmap(countsRaw, geneList, expFactors, plotDefs)
## pCountsHeatmapGrouped <- plotCountsHeatmapGrouped(countsRL, geneList, expFactors, plotDefs)
## plotCountsHeatmapGrouped(countsRL, geneList, expFactors, plotDefs)

# Specific analyses
heatmapList <- goHeatmaps(countsVSAvg, zscore = FALSE)
saveHeatmaps()

# Predicted activity or names for genes
# PF2 alignments to PTRM4 genome (also do with P. nodorum pf2)
# - Look at ones with more expression, more likely to be the actual transcription factor
# - Look at transcription factors
# - Map out more transcription factors for early and late (also fix gene clustering order, i.e. cluster at early but keep clustering order for late)
# - RNAi, chromatin remodelling (look at remodelling genes in other fungi + others), DNA methylation genes
# How are cultivars aligning? Magenta consistently with ZWB12, Wyalkatchem consistently with ZWB11
# Try with manual annotation of ORF for genes, or do with

# TODO: Heatmap with all SDEG, heatmap and SDEGs of between time factors

## n <- 10
## (countsRLAvg %>% scale())[
##   geneList$genes[names(geneList$genes) == c("ptrpf2", "pf2hyp1")], 1:3
## ]

# More differential expression: select groups

# Bar plots
for (t in 1:length(expFactors$factorList$time)) {
  counts <- t(scale(t(countsRL)))
  tstr <- expFactors$factorList$time[t]
  genes <- c(
    geneList$genes[["toxa"]],
    geneList$genes[["ptrpf2"]]
  )

  s <- c()
  for (i in 1:ncol(countsRL)) {
    if (length(grep(qq(".*@{toTitleCase(tstr)}.*"), expFactors$factorList$sample[i])) > 0) {
      s <- append(s, expFactors$factorList$sample[i])
    }
  }

  d <- data.frame()
  for (i in 1:length(genes)) {
    d <- rbind(d, data.frame(
      gene = genes[i],
      counts = counts[genes[i], s],
      group = rep(expFactors$factorList$cultivar, each = 3)
    ))
  }

  ggplot(d, aes(x = group, y = counts, fill = gene)) +
    geom_boxplot() +
    theme_std()
}

## pCountsHeatmapAll <- plotCountsHeatmap(countsRaw, geneList, expFactors, plotDefs, subset = FALSE, threshold = 1000)
pCountsHeatmapAll50 <- plotCountsHeatmapAll(countsRaw, expFactors, plotDefs, threshold = 50)

# Separate graphs across different factors
# 1: ONLY COMPARING AT EARLY TIME POINT ACROSS GROUPS
# 2: COMPARING BOTH TIME POINTS WITHIN GROUPS
p <- plotDiffExprTime(sdegGroups, sdeg)
## png(filename = "diffExprTime.png", width = 2500, height = 3500)
png(filename = "diffExprTime.png", width = 1200, height = 1800)
grid.arrange(p[[1]], p[[2]], p[[3]], heights = c(2, 2, 1))
dev.off()

# Information about SDEGs
# 11376 genes total, 6980 annotated
# 1099 SDEGs annotated of 1948 SDEGs
plotDiffExprAnnots(sdegGroups, sdeg$names)

# New clustered heatmap
## plotCountsHeatmapDiffExpr(countsRL, sdeg, plotDefs)

p1 <- grid.arrange(pTotalCounts, pTotalCountsBars, nrow = 1)
p2 <- grid.arrange(pTotalCountsNorm, pTotalCountsBarsNorm, nrow = 1)

savePlotGrid(p1, formatPath(expData[["species"]], "pTotalCounts_BOTH.png", root = FALSE), plotDefs)
savePlotGrid(p2, formatPath(expData[["species"]], "pTotalCountsNorm_BOTH.png", root = FALSE), plotDefs)
saveHeatmap(pCountsHeatmapAll50, formatPath(expData[["species"]], "pCountsHeatmapAll_top50.png", root = FALSE), plotDefs, width = 1000, height = 1000)
savePlotStd(pCountsHeatmap, formatPath(expData[["species"]], "pCountsHeatmap.png", root = FALSE), plotDefs, width = 6000, height = 4000)

# Specific genes
# Correlation: library size against actin expression
hkgCounts <- countsRaw[geneList[["genes"]][names(geneList[["genes"]]) == hkg], ]
pPairCountsScatter <- plotPairCountsScatter(hkgCounts, totalCountsAll, expFactors, plotDefs)
hkgCountsNorm <- countsRaw[geneList[["genes"]][names(geneList[["genes"]]) == hkg], ]
pPairCountsScatterNorm <- plotPairCountsScatter(hkgCountsNorm, totalCountsAllNorm, expFactors, plotDefs)

p3 <- grid.arrange(pPairCountsScatter, pPairCountsScatterNorm, nrow = 1)
savePlotGrid(p3, formatPath(expData[["species"]], "pPairCountsScatter_BOTH.png", root = FALSE), plotDefs, width = 6000, height = 4000)

savePlotGrid(
  ggarrange(pTotalCounts, pTotalCountsNorm, pTotalCountsBarsNorm, pPairCountsScatter, labels = c("A", "B", "C", "D")),
  formatPath(expData[["species"]], "paper1.png", root = FALSE),
  plotDefs,
  width = 6000,
  height = 6000,
)

# Alignment statistics
fileStarPtr <- formatPath(conf[["data"]][["data"]], conf[["data"]][["star"]], "ptr.tsv")
fileStarWheat <- formatPath(conf[["data"]][["data"]], conf[["data"]][["star"]], "wheat.tsv")

group <- rep(expFactors[["factorList"]][["group"]], each = 2)
for (i in 1:length(group)) {
  if (i %% 2 == 0) {
    group[i] <- paste0(group[i], "Disease")
  } else {
    group[i] <- paste0(group[i], "Control")
  }
}
groupOrdered <- rep(expFactors[["factorList"]][["groupOrdered"]], each = 2)
for (i in 1:length(groupOrdered)) {
  if (i %% 2 == 0) {
    groupOrdered[i] <- paste0(groupOrdered[i], "Disease")
  } else {
    groupOrdered[i] <- paste0(groupOrdered[i], "Control")
  }
}

## samplesPtr <- gsub("[^0-9]", "", expFactors[["factorList"]][["sample"]])
## samplesPtr <- gsub("1156", "", samplesPtr)
## samplesPtr <- gsub("1289", "", samplesPtr)
## samplesPtr <- as.numeric(samplesPtr)

starPtr <- read.table(fileStarPtr, header = TRUE)
starPtr <- starPtr[-(nrow(starPtr)), ]
starPtr[, "sample"] <- as.factor(starPtr[, "sample"])
starPtr[, "group"] <- rep(group, each = nreps)
starPtr[, "group"] <- factor(starPtr[, "group"], levels = groupOrdered)

cat <- colnames(starPtr)[c(-1, -length(starPtr))]
for (i in 1:nrow(starPtr)) {
  starPtr[i, cat] <- starPtr[i, cat] / starPtr[i, "input_reads"]
  starPtr[i, cat] <- 100 * starPtr[i, cat]
}
starPtrReads <- reshape2::melt(
  starPtr[, c("group", "mapped_unique", "mapped_multiple", "mapped_multiple_excess")]
)

pStarPtr <- ggplot(starPtrReads, aes(x = group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_std(plotDefs) +
  theme(
    axis.text.x = element_text(hjust = 1.0, vjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_discrete(labels = c("Uniquely mapping", "Multiple-mapping", "Excessively-mapping")) +
  ggtitle("Reads aligned to Ptr reference genome") +
  labs(x = "Group", y = "Percent (%)") +
  ylim(c(0, 80)) +
  guides(fill = guide_legend("Reads"))

starWheat <- read.table(fileStarWheat, header = TRUE)
starWheat <- starWheat[-(nrow(starWheat)), ]
starWheat[, "sample"] <- as.factor(starWheat[, "sample"])
starWheat[, "group"] <- rep(group, each = nreps)
starWheat[, "group"] <- factor(starWheat[, "group"], levels = groupOrdered)
cat <- colnames(starWheat)[c(-1, -length(starWheat))]

for (i in 1:nrow(starWheat)) {
  starWheat[i, cat] <- starWheat[i, cat] / starWheat[i, "input_reads"]
  starWheat[i, cat] <- 100 * starWheat[i, cat]
}
starWheatReads <- reshape2::melt(
  starWheat[, c("group", "mapped_unique", "mapped_multiple", "mapped_multiple_excess")]
)

pStarWheat <- ggplot(starWheatReads, aes(x = group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_std(plotDefs) +
  theme(
    axis.text.x = element_text(hjust = 1.0, vjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_discrete(labels = c("Uniquely mapping", "Multiple-mapping", "Excessively-mapping")) +
  ggtitle("Reads aligned to reference wheat reference genome") +
  labs(x = "Sample", y = "Percent (%)") +
  ylim(c(0, 80)) +
  guides(fill = guide_legend("Reads"))

pStarLegend <- cowplot::get_legend(pStarPtr + theme(legend.position = "right"))
pStarAlignment <- grid.arrange(
  ggarrange(pStarWheat, pStarPtr, labels = c("A", "B")),
  right = pStarLegend
)
savePlotGrid(pStarAlignment, formatPath(species, "paper2.png", root = FALSE), plotDefs)

# DEPRECATED: See above for new differential analysis
## # Venn diagrams
## results <- getResultsList(dds, expData, expFactors)
## sdegList <- getDifferentialGenes(results, expFactors)
## plotListAll <- plotDifferentialGenes(sdegList, expFactors, plotDefs, offsetR = 0, offsetC = 0, debug = FALSE)
## plotList <- plotListAll[[1]]
## legendList <- plotListAll[[2]]
## saveVennDiagrams(species, plotList, legendList, expFactors, plotDefs)


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
for (time in 1:length(expFactors[["factorList"]][["time"]])) {
  timeIndex <- time
  timePoint <- expFactors[["factorList"]][["time"]][time]

  modules <- getModules(clustRuns)
  modules <- modules[[timePoint]]
  goMap <- readMappings(fileGoMap)

  # Map cluster genes to their annotations using topGO
  # Reversed order of iterators... need to reinitialise lists
  clustPlots <- list()
  clustPlotsTG <- list()

  clustPlotsTG <- plotExprProfilesGo(clustPlotsTG, goMap, globalColours = TRUE)
  clustPlots <- plotExprProfilesGoRestricted(clustPlots, goMap, globalColours = TRUE, "aspergillus")

  # Aspergillus
  # Line graphs for each cluster
  # 1. Get subset of genes for each cluster
  # 2. Calculate std. dev of expression of genes
  # 3. Read in eigengene values for each cluster and cultivar
  # 4. Calculate z-score
  clustData <- read.table(
    formatPath(conf[["data"]][["data"]], conf[["data"]][["clust"]], species, timePoint, "Processed_Data", paste0(species, "RL.tsv_processed.tsv")),
    header = TRUE
  )
  rownames(clustData) <- clustData[, 1]
  colnames(clustData)
  clustData <- clustData[, -1]

  groups <- expFactors[["factorList"]][["group"]]
  groups <- groups[2 * 1:(length(groups) / 2) + timeIndex - 2]

  groupsOrdered <- expFactors[["factorList"]][["groupOrdered"]]
  groupsOrdered <- groupsOrdered[2 * 1:(length(groupsOrdered) / 2) + timeIndex - 2]

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
      theme_std_highres(plotDefsObj) +
      theme(
        axis.title.x = element_blank(),
      ) +
      # xlab(plotDefs[["defs"]][["names.xlab.cultivar"]]) +
      # ylab(plotDefs[["defs"]][["names.ylab.zscore"]]) +
      ggtitle(paste("Module", i))

    plotListME[[i]] <- p
  }

  ## plotListComplex1 <- list()
  ## for (i in 1:length(modules)) {
  ##   plotListComplex1[[i]] <- list()
  ##   for (j in 1:length(ont)) {
  ##     p <- arrangeGrob(clustPlotsTG[[i]][[j]][[2]], clustPlotsTG[[i]][[j]][[1]], plotListME[[i]], ncol = 1, heights = c(0.75, 0.25, 5.0))
  ##     savePlotGrid(
  ##       p,
  ##       formatPath(species, "clust", tolower(ont[j]), paste0("clust_goDist_module", i - 1, ".png"), root = FALSE),
  ##       plotDefs,
  ##       width = 6000,
  ##       height = 4000,
  ##     )
  ##     plotListComplex1[[i]][[j]] <- p
  ##   }
  ## }

  plotListComplex1 <- list()
  for (i in 1:length(modules)) {
    p <- ggarrange(
      ggarrange(
        grid.arrange(clustPlots[[i]][[1]][[2]], clustPlots[[i]][[1]][[1]], nrow = 2, heights = c(0.75, 0.25)),
        grid.arrange(clustPlots[[i]][[2]][[2]], clustPlots[[i]][[2]][[1]], nrow = 2, heights = c(0.75, 0.25)),
        grid.arrange(clustPlots[[i]][[3]][[2]], clustPlots[[i]][[3]][[1]], nrow = 2, heights = c(0.75, 0.25)),
        nrow = 3
      ),
      plotListME[[i]],
      nrow = 2,
      heights = c(0.6, 0.4),
      labels = i
    )
    plotListComplex1[[i]] <- p
  }

  plotListComplex2 <- list()
  for (i in 1:length(modules)) {
    p <- ggarrange(
      ggarrange(
        grid.arrange(clustPlotsTG[[i]][[1]][[2]], clustPlotsTG[[i]][[1]][[1]], nrow = 2, heights = c(0.75, 0.25)),
        grid.arrange(clustPlotsTG[[i]][[2]][[2]], clustPlotsTG[[i]][[2]][[1]], nrow = 2, heights = c(0.75, 0.25)),
        grid.arrange(clustPlotsTG[[i]][[3]][[2]], clustPlotsTG[[i]][[3]][[1]], nrow = 2, heights = c(0.75, 0.25)),
        nrow = 3
      ),
      plotListME[[i]],
      nrow = 2,
      heights = c(0.6, 0.4),
      labels = i
    )
    plotListComplex2[[i]] <- p
  }

  savePlotGrid(
    do.call(ggarrange, c(plotListComplex1, nrow = 1)) + theme_minimal(),
    formatPath(species, "clust", timePoint, "clust_goDist_aspergillus_all.png", root = FALSE),
    plotDefs,
    width = 6000 * length(modules),
    height = 6000
  )

  savePlotGrid(
    do.call(ggarrange, c(plotListComplex2, nrow = 1)) + theme_minimal(),
    formatPath(species, "clust", timePoint, "clust_goDist_all.png", root = FALSE),
    plotDefs,
    width = 6000 * length(modules),
    height = 6000
  )
}

# Aspergillus annotations
plots <- list()
for (i in 1:length(modules)) {
  plots[[i]] <- plotListComplex1[[i]]
}

p <- do.call(
  arrangeGrob,
  c(
    plots,
    nrow = 1
    ## top = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.lines.clustExpr.title\"]]}, Aspergillus DB [@{ont[i]}]")),
    ## bottom = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.xlab.cultivar\"]]}"))
  )
)
savePlotGrid(
  p,
  formatPath(species, "clust", paste0("clust_goDist_aspergillus_all", ".png"), root = FALSE),
  plotDefs,
  height = dim(p)[1] * 6000,
  width = dim(p)[2] * 6000
)

plots <- list()
for (i in 1:length(modules)) {
  plots[[i]] <- plotListComplex2[[i]]
}

# All GO annotations
p <- do.call(
  arrangeGrob,
  c(
    plots,
    nrow = 1
    ## top = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.lines.clustExpr.title\"]]}, Aspergillus DB [@{ont[i]}]")),
    ## bottom = textGrob(qq("@{plotDefs[[\"defs\"]][[\"names.xlab.cultivar\"]]}"))
  )
)
savePlotGrid(
  p,
  formatPath(species, "clust", paste0("clust_goDist_all", ".png"), root = FALSE),
  plotDefs,
  height = dim(p)[1] * 6000,
  width = dim(p)[2] * 6000
)

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
