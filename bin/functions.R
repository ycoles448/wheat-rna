# Functions
# TODO: Refactor into separate functions
appendGene <- function(gene = gene, name = name, geneList = genes) {
  geneList <- append(gene, geneList) # Actin
  names(geneList)[1] <- name

  return(geneList)
}

getSDEG <- function(contrasts = contrasts) {
  for (i in 1:length(contrasts)) {
    x <- padj
    y <- l2fc

    contrasts[[i]] <- subset(contrasts[[i]], padj < x & abs(log2FoldChange) > y)
  }

  return(contrasts)
}

plotContrastsHeatmap <- function(contrasts = contrasts) {
  plots <- list()
  data <- list()
  for (i in 1:length(contrasts)) {
    a <- padj
    b <- l2fc

    up <- data.frame(
      expr = nrow(subset(contrasts[[i]], padj < a & log2FoldChange > +b)),
      x = factor("Up")
    )
    dn <- data.frame(
      expr = nrow(subset(contrasts[[i]], padj < a & log2FoldChange < -b)),
      x = factor("Dn")
    )

    data[[i]] <- rbind(up, dn)
  }

  n <- c()
  for (i in data) {
    n <- append(n, i[1, "expr"])
    n <- append(n, i[2, "expr"])
  }

  for (i in 1:length(contrasts)) {
    d <- data[[i]]
    name <- cultivars[i]

    k <- c("Up", "Down")
    d[, "x"] <- factor(k, levels = k[2:1])

    p <- ggplot(d, aes(x = x, y = 0, fill = expr, width = 1, height = 1)) +
      geom_tile() +
      geom_text(aes(label = expr), size = 8) +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none"
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0, 0)) +
      ggtitle(name) +
      xlab("Expression") +
      coord_fixed() +
      scale_fill_gradient(low = "white", high = "#eb6e34", limits = c(min(n), max(n)))

    plots[[i]] <- p
  }

  return(plots)
}

saveContrastsHeatmap <- function(plots = plots) {
  ggsave(
    filename = "heatmap-FUNGI.png",
    grid.arrange(grobs = plots, widths = rep(1, times = length(plots))),
    width = 1000 * length(plots), height = 1000,
    units = "px",
    limitsize = FALSE
  )
  dev.off()
}

plotContrastsVolcano <- function(contrasts = contrasts) {
  plots <- list()
  data <- list()
  for (i in 1:length(contrasts)) {
    a <- padj
    b <- l2fc

    up <- data.frame(
      expr = nrow(subset(contrasts[[i]], padj < a & log2FoldChange > +b)),
      x = factor("Up")
    )
    dn <- data.frame(
      expr = nrow(subset(contrasts[[i]], padj < a & log2FoldChange < -b)),
      x = factor("Dn")
    )

    data[[i]] <- rbind(up, dn)
  }

  n <- c()
  for (i in data) {
    n <- append(n, i[1, "expr"])
    n <- append(n, i[2, "expr"])
  }

  for (i in 1:length(contrasts)) {
    d <- data[[i]]
    name <- cultivars[i]

    k <- c("Up", "Down")
    d[, "x"] <- factor(k, levels = k[2:1])

    p <- ggplot(d, aes(x = x, y = 0, fill = expr, width = 1, height = 1)) +
      geom_tile() +
      geom_text(aes(label = expr), size = 8) +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none"
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0, 0)) +
      ggtitle(name) +
      xlab("Expression") +
      coord_fixed() +
      scale_fill_gradient(low = "white", high = "#eb6e34", limits = c(min(n), max(n)))

    plots[[i]] <- p
  }

  return(plots)
}

allCounts <- function(counts = counts) {
  x <- c()
  for (i in 1:ncol(counts)) {
    x[i] <- sum(counts[, i])
  }
  names(x) <- colnames(counts)

  return(x)
}

averageCounts <- function(counts = counts) {
  cols <- ncol(counts)
  rows <- nrow(counts[genes, ])
  avgCounts <- data.frame(matrix(nrow = rows, ncol = cols / 3))

  for (i in 3 * (1:(cols / 3))) {
    for (r in 1:rows) {
      rownames(avgCounts)[r] <- rownames(counts[genes, ])[r]
      s <- sum(counts[genes, c(i - 2, i - 1, i)][r, ]) / 3
      avgCounts[r, (i / 3)] <- s
    }
  }
  colnames(avgCounts) <- meta[3 * (1:(cols / 3)), CONTRAST]

  return(avgCounts)
}

plotGenesHeatmap <- function(counts = counts, totalCounts = totalCounts, mode = "sample") {
  mode <- sym(mode)
  d1 <- reshape2::melt(as.matrix(totalCounts))
  d2 <- reshape2::melt(as.matrix(counts[genes, ]))
  d1 <- cbind(rep(cultivars, each = 2 * replicates), rep(groups, each = replicates), d1[, c(1, 3)])
  d2 <- cbind(rep(cultivars, each = 2 * replicates * length(genes)), rep(groups, each = replicates * length(genes)), d2[, c(2, 3, 1)])
  d2 <- cbind(
    d2,
    rep(names(genes), (nrow(d2) / length(genes)), length(genes))
  )
  colnames(d1) <- c("cultivar", "group", "sample", "count")
  colnames(d2) <- c("cultivar", "group", "sample", "count", "gene", "cat")
  d1[, "sample"] <- samples
  d2[, "sample"] <- rep(samples, each = length(genes))
  d1[, "group"] <- factor(d1[, "group"], levels = groupsOrdered)
  d2[, "group"] <- factor(d2[, "group"], levels = groupsOrdered)
  d1[, "sample"] <- factor(d1[, "sample"], levels = samplesOrdered)
  d2[, "sample"] <- factor(d2[, "sample"], levels = samplesOrdered)

  # TODO: Vary colours of genes by category
  # Assign discrete colours using HSV
  max <- max(d2[, "count"])
  min <- min(d2[, "count"])
  k <- max - min
  h <- c()
  for (i in 1:length(unique(d2[, "cat"]))) {
    h[i] <- 0 + i * (1 / length(unique(d2[, "cat"])))
  }
  names(h) <- unique(d2[, "cat"])
  for (i in 1:nrow(d2)) {
    ## message(paste("Row:", i))
    val <- d2[i, "count"]
    cat <- d2[i, "cat"]

    hue <- h[which(names(h) == cat)]
    sat <- sqrt(val / k)
    val <- 1.0

    ## print(paste0("Gene: ", d2[i, "gene"]))
    ## print(paste0("> Hue: ", hue))
    ## print(paste0("> Sat: ", sat))
    ## print(paste0("> Val: ", val))

    # The dumb conversion: RGB -> HSV for colour blending
    x <- rgb(1, 0, 0)
    r1 <- sqrt(strtoi(paste0("0x", substr(x, 2, 3))) / 255)
    g1 <- sqrt(strtoi(paste0("0x", substr(x, 4, 5))) / 255)
    b1 <- sqrt(strtoi(paste0("0x", substr(x, 6, 7))) / 255)

    y <- rgb2hsv(r = r1, g = g1, b = b1)
    ratio <- 10
    clr <- hsv(((ratio - 1) * y[1] + hue) / ratio, sat, val)

    d2[i, "colour"] <- clr
  }

  margins <- unit(c(0, 1, 0, 1), "cm")
  margins.axisY <- element_text(margin = unit(c(0, 0, 0, 1), "cm"))

  p1 <- ggplot(d1, aes(x = !!mode, y = count, fill = cultivar)) +
    theme_classic() +
    theme(
      ## axis.text.x = element_text(angle = 90, hjust = 1)
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = margins.axisY,
      axis.ticks.length.x = unit(0, "cm"),
      plot.margin = margins,
      legend.position = "right"
    ) +
    guides(fill = "none") +
    geom_col() +
    ggtitle(paste0("Raw read counts per ", mode))

  p2 <- ggplot(d2, aes(x = !!mode, y = gene, fill = colour)) +
    geom_tile() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = margins.axisY,
      plot.margin = margins,
      legend.position = "right"
    ) +
    ## guides(colour = "legend") +
    guides(fill = "none") +
    scale_fill_identity() +
    ggtitle(paste0("Gene counts across ", mode, "s"))
  p2

  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g <- rbind(g1, g2, size = "first")
  grid.arrange(g)

  ggsave(
    filename = "reads+expression+sqrt.png",
    grid.arrange(g),
    width = 1000, height = 2000,
    units = "px",
    limitsize = FALSE
  )
  dev.off()

  return(list(plot = g, d1 = d1, d2 = d2, p1 = p1, p2 = p2))
}

plotGeneBox <- function(gene = NA, counts = counts) {
  geneName <- toTitleCase(toTitleCase(names(gene)[1]))
  d <- reshape2::melt(counts[gene, ])
  d[, "cultivar"] <- rep(cultivars, each = 2 * replicates)
  d[, "sample"] <- samples
  d[, "group"] <- rep(groups, each = replicates)
  d[, "group"] <- factor(d[, "group"], levels = groupsOrdered)
  colnames(d)[1] <- "count"

  ## d <- d[, c("group", "sample", "count")]
  p <- ggplot(d, aes(x = group, y = count, fill = cultivar)) +
    geom_boxplot() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    guides(fill = "none") +
    ggtitle(paste0("Expression: ", geneName))

  return(p)
}

plotGeneNormBox <- function(gene = NA, counts = counts, norm = NA) {
  geneName <- toTitleCase(names(gene)[1])
  normName <- toTitleCase(names(norm)[1])
  d <- reshape2::melt(counts[gene, ] / counts[norm, ])
  d[, "cultivar"] <- rep(cultivars, each = 2 * replicates)
  d[, "sample"] <- samples
  d[, "group"] <- rep(groups, each = replicates)
  d[, "group"] <- factor(d[, "group"], levels = groupsOrdered)
  colnames(d)[1] <- "ratio"

  ## d <- d[, c("group", "sample", "count")]
  p <- ggplot(d, aes(x = group, y = ratio, fill = cultivar)) +
    geom_boxplot() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    guides(fill = "none") +
    ggtitle(paste0("Expression: ", geneName, " against ", normName))

  return(p)
}

getPairs <- function(dds = dds) {
  x <- max(as.numeric(sapply(meta[, CONTRAST], function(x) {
    substr(x, 1, 1)
  })))

  pairs <- list()
  k <- 1
  for (time in times) {
    c <- 1
    for (i in 1:x) {
      for (j in 2:x) {
        if (j > c) {
          # Indexes are (i, j)
          # Pairwise comparisons of disease against controls
          g1 <- paste0(i, toTitleCase(controls[2]), toTitleCase(time))
          g2 <- paste0(j, toTitleCase(controls[1]), toTitleCase(time))

          dds.res <- results(
            dds,
            alpha = ALPHA,
            contrast = c(CONTRAST, g1, g2),
            filterFun = ihw,
            parallel = TRUE
          )
          pairs[[k]] <- dds.res
          k <- k + 1
        }
      }
      c <- c + 1
    }
  }
  return(pairs)
}

getGenes <- function(pairs = pairs) {
  genes <- list()
  for (i in 1:length(pairs)) {
    dds <- pairs[[i]]
    up <- subset(dds, log2FoldChange > +L2FC & padj < ALPHA)
    dn <- subset(dds, log2FoldChange < -L2FC & padj < ALPHA)
    genes[[i]] <- list(up = up, dn = dn)
  }

  return(genes)
}

getCommonGenes <- function(genes = genes) {
  # Get upregulated genes across each cultivar (common to cultivars)
  # Late infection
  allUp <- list()
  allDn <- list()
  c <- 1
  for (i in (length(genes) / 2 + 1):length(genes)) {
    up <- rownames(genes[[i]][["up"]])
    dn <- rownames(genes[[i]][["dn"]])
    allUp[[c]] <- up
    allDn[[c]] <- dn
    c <- c + 1
  }

  return(list(up = Reduce(intersect, allUp), dn = Reduce(intersect, allDn)))
}

writeRegGenes <- function(x = OUT, genes = commonGenes) {
  if (!dir.exists(x)) {
    dir.create(x)
  }

  for (i in c("up", "dn")) {
    fwrite(x = as.list(genes[[i]]), file = paste(x, paste0(i, "SDEG.tsv"), sep = "/"), sep = "\n")
  }
}

analyseGO <- function(map = MAP, genes = genes) {
  # GO analysis is per group, use common significant genes per group
  goMap <- readMappings(MAP)
  geneNames <- factor(row.names(genes[[1]][["up"]]))

  sigGenes <- factor(as.integer(names(goMap) %in% geneNames))
  names(sigGenes) <- names(goMap)

  # Expression matrix
  x <- genes[[1]][["up"]][, "padj"]
  x <- -log10(x)
  x <- x / max(x)
  names(x) <- rownames(genes[[1]][["up"]])

  # Empty matrix of all genes
  y <- rep(0, times = length(goMap))
  names(y) <- names(goMap)

  # Get genes in expression matrix only in empty matrix
  z <- x[names(x) %in% names(goMap)]
  for (gene in names(z)) {
    y[gene] <- z[gene]
  }

  go <- new(
    "topGOdata",
    ontology = "MF",
    allGenes = sigGenes,
    geneSel = factor(rownames(genes[[1]][["up"]])),
    annot = annFUN.gene2GO,
    gene2GO = goMap
  )

  # Iterator
  for (k in c(0, length(pairs) / 2)) {
    x <- length(pairs) / 4
    j <- 1
    for (i in (1 + k):(length(pairs) / 2 + k)) {
      print(paste("i:", i, "x:", x, "j:", j))
      for (e in c("up", "dn")) {

      }

      if (j != x) {
        j <- j + 1
      } else {
        x <- x - 1
        j <- 1
      }
    }
  }
}
