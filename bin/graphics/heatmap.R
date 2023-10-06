plotCountsHeatmap <- function(counts, geneListObj, expObjFactors, plotDefsObj, subset = TRUE) {
  if (subset == TRUE) {
    counts <- counts[geneListObj[["genes"]], ]
  }

  d <- reshape2::melt(counts)
  defs <- plotDefsObj[["defs"]]
  colnames(d) <- c("gene", "sample", "count")

  countMax <- max(counts)
  countMin <- min(counts)

  # Assign colours to each gene, based on count and group
  d[, "cat"] <- names(geneListObj[["genes"]])
  d[, "colour"] <- "#000000"
  d[, "group"] <- factor(
    rep(
      expObjFactors[["factorList"]][["group"]],
      each = length(geneListObj[["genes"]])
    ),
    levels = expObjFactors[["factorList"]][["groupOrdered"]]
  )
  range <- countMax - countMin
  cats <- unique(names(geneListObj[["genes"]]))
  hues <- vector(length = length(cats))
  names(hues) <- cats
  ratio <- 5
  for (i in 1:length(cats)) {
    hues[i] <- i * (1 / length(cats))
  }

  for (i in 1:nrow(d)) {
    h1 <- hues[d[, "cat"][i]]
    c2 <- clr2rgb(defs["clr2"])
    h2 <- rgb2hsv(c2[1], c2[2], c2[3])[1]

    h <- ((ratio - 1) * h2 + h1) / ratio
    s <- (d[i, "count"] / range)^(1 / 3)
    v <- 1.0

    d[i, "colour"] <- hsv(h, s, v)
  }

  breaks <- c()
  for (i in 1:length(cats)) {
    # Get maximum value for each category
    h <- ((ratio - 1) * h2 + hues[i]) / ratio
    s <- (max(d[d[, "cat"] == cats[i], "count"]) / range)^(1 / 3)
    v <- 1.0

    breaks[i] <- hsv(h, s, v)
  }

  d[, "gene"] <- factor(
    geneList[["genes"]],
    levels = geneList[["genes"]]
  )

  ## print(breaks)
  ## print(cats)

  p <- ggplot(d, aes(x = group, y = gene, fill = colour)) +
    geom_tile() +
    ggtitle(defs["names.heatmap.geneCounts.title"]) +
    xlab(defs["names.xlab.cultivar"]) +
    ylab(defs["names.ylab.totalCounts"]) +
    theme_std(plotDefsObj) +
    guides(
      group = expObjFactors[["factorList"]][["group"]]
    ) +
    scale_fill_identity(
      defs["names.heatmap.geneCounts.legend"],
      labels = factor(cats, levels = rev(cats)),
      breaks = factor(breaks),
      guide = "legend"
    )

  return(p)
}

plotCountsHeatmapAll <- function(counts, expObjFactors, plotDefsObj, threshold = 100) {
  x <- vector()
  for (row in 1:nrow(counts)) {
    x[row] <- sum(counts[row, ])
    names(x)[row] <- rownames(counts)[row]
  }

  colnames(counts) <- expObjFactors[["factorList"]][["sample"]]

  x <- names(tail(sort(x), n = threshold))
  counts <- counts[x, ]

  p <- Heatmap(
    counts,
    row_labels = x,
    use_raster = TRUE,
    raster_quality = 1,
    raster_device = "png",
  )
  return(p)
}
