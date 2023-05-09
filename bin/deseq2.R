# =========================== #
# - DESeq2 RNA-seq analysis - #
# =========================== #

# Requires matrices for gene counts and experimental design
# Author: Yutathkarn Coles, Paula Moolhuijzen


# Load packages into R environment
library("DESeq2")
library("BiocManager")
library("BiocParallel")
library("Rgraphviz")
library("IHW")
library("RColorBrewer")
library("ggplot2")
library("ggVennDiagram")
library("ggpubr")
library("topGO")
library("data.table")
library("stringr")
library("ballgown")


# Allow parallelisation
register(MulticoreParam(12))


# Plotting functions
# Emulation of default ggplot colours
ggColours <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

counts2plot <- function(x = NA, colour = NA) {
  p <- barplot(
    colSums(x) / 1000000,
    horiz = TRUE,
    col = colour,
    las = 1,
    xlab = "No. reads aligned (million)",
    cex.axis = 20,
    cex.names = 20,
    cex.lab = 20,
    border = NA
  )
  return(p)
}

# PCA
pca2plot <- function(output, pca = pca, title = NA, type = NA) {
  if (type == "fungi") {
    a <- aes(
      x = PC1, y = PC2,
      color = cultivar,
      label = name,
      shape = time,
      size = rating
    )
  } else if (type == "wheat") {
    a <- aes(
      x = PC1, y = PC2,
      color = cultivar,
      label = name,
      shape = control,
      size = rating
    )
  }

  p <- ggplot(
    pca,
  ) +
    a +
    scale_size_continuous(
      labels = c("R", "MR", "MS", "SVS"),
      range = c(5, 10)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    guides(shape = guide_legend(override.aes = list(size = 5))) +
    geom_point() +
    xlab(paste0("PC1: ", pct[1], "% var")) +
    ylab(paste0("PC2: ", pct[2], "% var")) +
    ggtitle(paste("PCA, variance stabilised", title, sep = ", ")) +
    theme_classic()
  return(p)
  ## ggsave(output, plot = p, units = "px", width = 3840, height = 2160)
}


# Loading data and adding metadata not present in metadata matrix
# TODO: Add data into matrix instead of adding data here
# Loading in data
## counts.wheat <- read.csv("wheat-matrix.tsv", sep = "\t", header = TRUE)
counts.wheat <- read.csv("wheat-matrix.unstranded.tsv", sep = "\t", header = TRUE)
counts.fungi <- read.csv("ptr-matrix.tsv", sep = "\t", header = TRUE)
meta.wheat <- read.csv("meta-reduced.tsv", sep = "\t", header = TRUE)

# Other useful data
# Cultivars are ordered by disease rating, and then alphabetically
cultivars <- c("Machete", "Phantom", "Yitpi", "Wyalkatchem", "Magenta", "ZWB11Q56", "ZWB12Q89")

# Adding metadata to samples
meta.wheat[, "toxa"] <- "no"
meta.wheat[, "toxc"] <- "unknown"
meta.wheat[, "rating"] <- "SVS"

# ToxA
meta.wheat[, "toxa"][meta.wheat[, "cultivar"] == "Yitpi"] <- "susceptible"
meta.wheat[, "toxa"][meta.wheat[, "cultivar"] == "Phantom"] <- "susceptible"

# ToxC
meta.wheat[, "toxc"][meta.wheat[, "cultivar"] == "ZWB12Q89"] <- "resistant"
meta.wheat[, "toxc"][meta.wheat[, "cultivar"] == "ZWB11Q56"] <- "resistant"
meta.wheat[, "toxc"][meta.wheat[, "cultivar"] == "Magenta"] <- "resistant"
meta.wheat[, "toxc"][meta.wheat[, "cultivar"] == "Phantom"] <- "susceptible"

# Rating
meta.wheat[, "rating"][meta.wheat[, "cultivar"] == "ZWB12Q89"] <- "R"
meta.wheat[, "rating"][meta.wheat[, "cultivar"] == "ZWB11Q56"] <- "R"
meta.wheat[, "rating"][meta.wheat[, "cultivar"] == "Magenta"] <- "MR"
meta.wheat[, "rating"][meta.wheat[, "cultivar"] == "Wyalkatchem"] <- "MS"

# Groups
c <- 1
for (i in cultivars) {
  meta.wheat[meta.wheat[, "cultivar"] == i & meta.wheat[, "time"] == "early", "group"] <- c
  c <- c + 1
  meta.wheat[meta.wheat[, "cultivar"] == i & meta.wheat[, "time"] == "late", "group"] <- c
  c <- c + 1
}
meta.wheat[, "group"] <- as.roman(meta.wheat[, "group"])


# Colours
c <- ggColours(7)
for (i in seq(1, 7)) {
  meta.wheat[meta.wheat[, "cultivar"] == cultivars[i], "colour"] <- c[i]
}

# Correcting columns to factors
for (i in c("sample", "cultivar", "control", "time", "treatment", "pool", "row", "toxa", "toxc", "rating", "group")) {
  meta.wheat[, i] <- as.factor(meta.wheat[, i])
}

# Remove controls from samples and metadata for fungi analysis
meta.fungi <- meta.wheat[meta.wheat[, "control"] != "control", ]
meta.fungi <- meta.fungi[, names(meta.fungi) != "control"]
counts.fungi <- counts.fungi[, colnames(counts.wheat) %in% meta.fungi[, "sample"]]


# Reads per sample


# Setting variables
# Set column to use for grouping
f <- "rating"
eff <- 1
pval <- 0.1

# Set contrast
# First two contrasts are compared in Venn diagram
contrasts <- list(
  c("rating", "R", "SVS"),
  c("time", "early", "late"),
  c("toxc", "resistant", "susceptible")
)


## # Fungi analysis
## dds.fungi <- DESeqDataSetFromMatrix(
##   countData = counts.fungi,
##   colData = meta.fungi,
##   design = formula("~ time + cultivar")
## )
## dds.fungi <- DESeq(dds.fungi, parallel = TRUE)
## counts.fungi.norm <- counts(dds.fungi, normalized = TRUE)

## vst.fungi <- vst(dds.fungi, blind = FALSE)
## pca.fungi <- plotPCA(vst.fungi, intgroup = f, returnData = TRUE, ntop = 2000)
## pct <- round(100 * attr(pca.fungi, "percentVar"))

## pca.fungi[, "cultivar"] <- meta.fungi[, "cultivar"]
## pca.fungi[, "rating"] <- meta.fungi[, "rating"]
## pca.fungi[, "rating"] <- rep(c(4, 2, 4, 3, 4, 1, 1), each = 6)
## pca.fungi[, "time"] <- meta.fungi[, "time"]

## dds.fungi <- DESeq(dds.fungi, parallel = TRUE)
## res.fungi <- results(
##   dds.fungi,
##   alpha = 0.05,
##   contrast = unlist(contrasts[1]),
##   filterFun = ihw,
##   parallel = TRUE
## )
## sig.fungi <- subset(res.fungi, (padj < 0.05 & !is.na(padj)) & abs(log2FoldChange) >= eff)


# Wheat analysis
# Can compare wheat samples by:
#   1. Disease rating
#   2. Cultivar
#   3. Controls
dds.wheat <- DESeqDataSetFromMatrix(
  countData = counts.wheat,
  colData = meta.wheat,
  design = formula("~ time + control + rating")
)

vst.wheat <- vst(dds.wheat, blind = FALSE)
pca.wheat <- plotPCA(vst.wheat, intgroup = f, returnData = TRUE, ntop = 2000)
pct <- round(100 * attr(pca.wheat, "percentVar"))

# Add metadata to PCA, factors used for ggplot
pca.wheat[, "control"] <- meta.wheat[, "control"]
pca.wheat[, "cultivar"] <- meta.wheat[, "cultivar"]
pca.wheat[, "rating"] <- rep(c(4, 2, 4, 3, 4, 1, 1), each = 12)
pca.wheat[, "time"] <- meta.wheat[, "time"]

## p1 <- pca2plot(output = "fungi-pca.png", title = "fungi", pca = pca.fungi, type = "fungi")
## p2 <- pca2plot(output = "wheat-pca.png", title = "wheat", pca = pca.wheat, type = "fungi")
## p3 <- pca2plot(output = NA, pca = pca.wheat[pca.wheat[, "control"] == "control", ], title = "wheat controls", type = "fungi")
## p4 <- pca2plot(output = NA, pca = pca.wheat[pca.wheat[, "control"] == "disease", ], title = "wheat disease", type = "fungi")

## ggsave(filename = "wheat+fungi-pca.png", plot = ggarrange(p1, p2), width = 6000, height = (6000 / 16) * 9, units = "px")
## ggsave(filename = "wheat-controls+disease-pca.png", plot = ggarrange(p3, p4), width = 6000, height = (6000 / 16) * 9, units = "px")

dds.wheat <- DESeq(dds.wheat, parallel = TRUE)
res.wheat <- results(
  dds.wheat,
  alpha = 0.05,
  contrast = unlist(contrasts[1]),
  filterFun = ihw,
  parallel = TRUE
)
sig.wheat <- subset(res.wheat, (padj < pval & !is.na(padj)) & abs(log2FoldChange) >= eff)


## sig.up <- subset(sig, log2FoldChange >= 1)
## sig.dn <- subset(sig, log2FoldChange <= 1)


## for (c in contrasts) {
##   dds <- DESeqDataSetFromMatrix(
##     countData = counts,
##     colData = meta,
##     ## design = ~ cultivar + time + pool + row + toxa + toxc + control
##     design = formula(paste0("~ time + ", c[1], " + control"))
##   )

##   # Filtering of counts
##   # Normalisation of counts
##   dds <- DESeq(dds, parallel = TRUE)
##   counts.norm <- counts(dds, normalized = TRUE)

##   # Filtering out reads with low counts
##   keep <- rowSums(counts(dds)) >= 0
##   dds <- dds[keep, ]
##   res <- results(
##     dds,
##     alpha = 0.05,
##     contrast = c,
##     filterFun = ihw,
##     parallel = TRUE
##   )

##   eff <- 1
##   sig <- subset(res, (padj <= 0.05 & !is.na(padj) & abs(log2FoldChange) >= eff))
##   sig.up <- subset(sig, log2FoldChange >= eff)
##   sig.dn <- subset(sig, log2FoldChange <= eff)

##   assign(paste("sig", c, sep = "."), sig)
##   assign(paste("sig.up", c, sep = "."), sig.up)
##   assign(paste("sig.dn", c, sep = "."), sig.dn)

##   # Plotting differential expression
##   ## plotMA(res, ylim = c(-15, 15))
##   ## plotMA(res.lfc, ylim = c(-15, 15))

##   ## plotcounts(
##   ##   dds,
##   ##   gene = which.min(res.lfc$padj),
##   ##   intgroup = "control")

##   ## n <- paste0(c[1], "-", c[2], "-", c[3], "-sig.tsv")
##   ## write.table(sig, n, sep = "\t")
## }


## # Compare significant genes
## l1 <- read.csv(
##   paste0(
##     unlist(contrasts[1])[1], "-",
##     unlist(contrasts[1])[2], "-",
##     unlist(contrasts[1])[3], "-sig.tsv"
##   ),
##   sep = "\t"
## )
## l2 <- read.csv(
##   paste0(
##     unlist(contrasts[2])[1], "-",
##     unlist(contrasts[2])[2], "-",
##     unlist(contrasts[2])[3], "-sig.tsv"
##   ),
##   sep = "\t"
## )

## l1.genes <- rownames(l1)
## l2.genes <- rownames(l2)

## venn <- list(l1.genes, l2.genes)
## p <- ggVennDiagram(
##   venn,
##   category.names = c(unlist(contrasts[1])[1], unlist(contrasts[2])[1])
## ) +
##   ggtitle(
##     paste0(
##       "Significant genes from contrasts\n",
##       unlist(contrasts[1])[1], " and ",
##       unlist(contrasts[2])[1]
##     )
##   )
## ## p


# GO term enrichment analysis
# g2go <- readMappings("data/ptr-m4.g2go.map")
g2go <- readMappings("data/wheat-cs.g2go.map")
genes <- names(g2go)

# List of genes for topGO are padj values named with gene names
terms <- c("CC", "BP", "MF")

sig.up <- subset(sig.wheat, log2FoldChange > eff)
sig.dn <- subset(sig.wheat, log2FoldChange < eff)

sig.up
sig.dn

l.up <- sig.up[, "pvalue"]
l.dn <- sig.dn[, "pvalue"]
names(l.up) <- rownames(sig.up)
names(l.dn) <- rownames(sig.dn)

selGenes <- function(x) {
  return(x < 0.01)
}

l <- list()
i <- 1
for (t in terms) {
  go.up <- new(
    "topGOdata",
    ontology = t,
    allGenes = l.up,
    geneSel = selGenes,
    annot = annFUN.gene2GO,
    gene2GO = g2go
  )

  go.dn <- new(
    "topGOdata",
    ontology = t,
    allGenes = l.dn,
    geneSel = selGenes,
    annot = annFUN.gene2GO,
    gene2GO = g2go
  )

  results.up.fish <- runTest(
    go.up,
    algorithm = "weight01",
    statistic = "fisher"
  )

  results.dn.fish <- runTest(
    go.dn,
    algorithm = "weight01",
    statistic = "fisher"
  )

  reg.up <- GenTable(
    go.up,
    weightFisher = results.up.fish,
    orderBy = "weightFisher",
    topNodes = length(results.up.fish@score),
    numChar = 120
  )

  reg.dn <- GenTable(
    go.dn,
    weightFisher = results.dn.fish,
    orderBy = "weightFisher",
    topNodes = length(results.dn.fish@score),
    numChar = 120
  )

  ## reg.up[, "weightFisher"] <- gsub("<", "", as.character(reg.up[, "weightFisher"]))
  ## reg.up[, "weightFisher"] <- as.numeric(reg.up[, "weightFisher"])
  ## reg.up.p <- reg.up[which(reg.up[, "weightFisher"] < 0.01), ]
  reg.up.p <- reg.up
  reg.up.p[, "reg"] <- "up"
  assign(paste("reg.up", t, sep = "."), reg.up)

  ## reg.dn[, "weightFisher"] <- gsub("<", "", as.character(reg.dn[, "weightFisher"]))
  ## reg.dn[, "weightFisher"] <- as.numeric(reg.dn[, "weightFisher"])
  ## reg.dn.p <- reg.dn[which(reg.dn[, "weightFisher"] < 0.01), ]
  reg.dn.p <- reg.dn
  reg.dn.p[, "reg"] <- "dn"
  assign(paste("reg.dn", t, sep = "."), reg.dn)

  go <- rbind(reg.up.p, reg.dn.p)
  go[, "reg"] <- factor(
    go[, "reg"],
    levels = c("up", "dn"),
    labels = c("Up", "Down")
  )
  go[, "domain"] <- t
  go[, "domain"] <- factor(
    go[, "domain"],
    levels = t,
    labels = t
  )

  l[[i]] <- go
  i <- i + 1
}


# Plot GO terms: metabolic function (MF)
go <- rbindlist(l, use.names = TRUE)
go[, "species"] <- "Wheat"
mf <- subset(go, domain == "MF")
mf <- subset(mf, mf$Significant / mf$Annotated > 0.7 & mf$weightFisher < 0.3)

n <- length(unique(mf$Term))
k <- mf$Significant
set <- colorRampPalette(brewer.pal(12, "Paired"))(n)

mf$Term <- str_wrap(mf$Term, 30)

text_std <- element_text(size = 12)

ggplot(mf, aes(species)) +
  theme_classic() +
  theme(
    legend.position = "right"
  ) +
  ## theme(
  ##   axis.text = text_std,
  ##   axis.title = text_std,
  ##   legend.text = text_std
  ## ) +
  geom_bar(data = subset(mf, reg == "Up"), aes(y = Significant, fill = Term), stat = "identity", position = "stack") +
  geom_bar(data = subset(mf, reg == "Down"), aes(y = -Significant, fill = Term), stat = "identity", position = "stack") +
  xlab("Species") +
  ylab("Number of SDEG regulated") +
  guides(fill = guide_legend(nrow = n, title = "GO molecular function")) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = set)


## for (t in terms) {
##   go.1 <- new("topGOdata", ontology = t, allGenes = )
## }

# Analysis using Ballgown
# Fungi analysis
## bg <- ballgown()
