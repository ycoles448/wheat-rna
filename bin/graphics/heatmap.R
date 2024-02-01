plotCountsHeatmapDiffExpr <- function(counts, sdeg, plotDefsObj = plotDefs) {
  c = counts[sdeg$names, ]

  xdist <- hclust(dist(t(c), method = "euclidean"), method = "ward.D")
  ydist <- hclust(dist(c, method = "euclidean"), method = "ward.D")

  d <- reshape2::melt(c)
  names(d) <- c("gene", "id", "count")
  d[, "id"] <- factor(as.character(d[, "id"]), levels = xdist$labels[xdist$order])
  d[, "gene"] <- factor(as.character(d[, "gene"]), levels = ydist$labels[ydist$order])

  p <- ggplot(d, aes(x = id, y = gene, fill = count)) +
    geom_tile() +
    ## ggdendrogram(xdist) +
    theme_std() +
    theme(
      axis.text.y = element_blank(),
      legend.position = "NONE",
    )
  p

  ## p <- heatmaply(d)
    ## theme_std() +
    ## theme(
    ##   axis.text.x = element_blank()
    ## )

  ## ggarrange(p, ggdendrogram(xdist), ncol = 1)
}

plotCountsHeatmap <- function(counts, geneListObj, expObjFactors, plotDefsObj, subset = TRUE) {
  if (subset == TRUE) {
    counts <- counts[geneListObj[["genes"]], ]
  }

  d <- reshape2::melt(counts)
  defs <- plotDefsObj[["defs"]]

  colnames(d) <- c("gene", "sample", "count")

  countMax <- max(counts)
  countMin <- min(counts)

  # Generate hues
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
  ## hues <- brewer.pal(length(cats), "Accent")
  hues <- rainbow(length(cats))
  names(hues) <- cats

  # Colour options
  ratio <- 1
  k <- 4

  ## for (i in 1:length(cats)) {
  ##   hues[i] <- i * (1 / length(cats))
  ## }

  for (i in 1:nrow(d)) {
    ## h1 <- hues[d[, "cat"][i]]
    c1 <- clr2rgb(hues[d[, "cat"][i]])
    h1 <- rgb2hsv(c1[1], c1[2], c1[3])[1]

    c2 <- clr2rgb(defs["clr2"])
    h2 <- rgb2hsv(c2[1], c2[2], c2[3])[1]

    ## h <- ((ratio - 1) * h2 + h1) / ratio
    ## s <- (d[i, "count"] / range)^(1 / 4)
    ## v <- 1.0

    h <- h1
    s <- (d[i, "count"] / range)^(1 / k)
    v <- 1.0

    d[i, "colour"] <- hsv(h, s, v)
    ## print(d[i, "colour"])
  }

  breaks <- c()
  for (i in 1:length(cats)) {
    # Get maximum value for each category
    ## h <- ((ratio - 1) * h2 + hues[i]) / ratio
    c <- clr2rgb(hues[i])
    h <- rgb2hsv(c[1], c[2], c[3])[1]
    s <- (max(d[d[, "cat"] == cats[i], "count"]) / range)^(1 / k)
    v <- 1.0

    breaks[i] <- hsv(h, s, v)
  }

  d[, "gene"] <- factor(
    geneList[["genes"]],
    ## levels = geneList[["genes"]]
    levels = unique(geneList[["genes"]])
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

plotCountsHeatmapGrouped <- function(counts, geneListObj, expObjFactors, plotDefsObj, subset = TRUE) {
  if (subset == TRUE) {
    counts <- counts[geneListObj[["genes"]], ]
  }

  ## dclust <- hclust(dist(counts, method = "euclidean"), method = "ward.D")
  xdist <- hclust(dist(t(counts), method = "euclidean"), method = "ward.D")
  ydist <- hclust(dist(counts, method = "euclidean"), method = "ward.D")

  ## Only plotting the first 3 supplied categories?
  d <- reshape2::melt(counts)
  defs <- plotDefsObj[["defs"]]
  colnames(d) <- c("gene", "sample", "count")

  countMax <- max(counts)
  countMin <- min(counts)
  ## cats <- unique(names(geneListObj[["genes"]]))

  # Generate hues
  # Assign colours to each gene, based on count and group
  ## d[, "cat"] <- names(geneListObj[["genes"]])
  d[, "sample"] <- factor(as.character(d[, "sample"]), levels = xdist$labels[xdist$order])
  d[, "gene"] <- factor(as.character(d[, "gene"]))

  p <- ggplot(d, aes(x = sample, y = gene, fill = count)) +
    geom_tile() +
    ggtitle(defs["names.heatmap.geneCounts.title"]) +
    xlab(defs["names.xlab.cultivar"]) +
    ylab(defs["names.ylab.totalCounts"]) +
    theme_std(plotDefsObj) +
    theme(
      ## axis.text.y = element_blank(),
    ) +
    guides(
      group = expObjFactors[["factorList"]][["group"]]
    )
    ## scale_fill_identity(
    ##   defs["names.heatmap.geneCounts.legend"],
    ##   labels = factor(cats, levels = rev(cats)),
    ##   breaks = factor(breaks),
    ##   guide = "legend"
    ## )

  return(p)
}

plotCountsHeatmapGroupedSubset <- function(counts, geneVec, expObjFactors, plotDefsObj) {
  counts <- counts[geneVec, ]

  ## dclust <- hclust(dist(counts, method = "euclidean"), method = "ward.D")
  xdist <- hclust(dist(t(counts), method = "euclidean"), method = "ward.D")
  ydist <- hclust(dist(counts, method = "euclidean"), method = "ward.D")

  ## Only plotting the first 3 supplied categories?
  d <- reshape2::melt(counts)
  defs <- plotDefsObj[["defs"]]
  colnames(d) <- c("gene", "sample", "count")

  countMax <- max(counts)
  countMin <- min(counts)
  ## cats <- unique(names(geneListObj[["genes"]]))

  # Generate hues
  # Assign colours to each gene, based on count and group
  ## d[, "cat"] <- names(geneListObj[["genes"]])
  d[, "sample"] <- factor(as.character(d[, "sample"]), levels = xdist$labels[xdist$order])
  d[, "gene"] <- factor(as.character(d[, "gene"]), levels = ydist$labels[ydist$order])

  p <- ggplot(d, aes(x = sample, y = gene, fill = count)) +
    geom_tile() +
    ggtitle(defs["names.heatmap.geneCounts.title"]) +
    xlab(defs["names.xlab.cultivar"]) +
    ylab(defs["names.ylab.totalCounts"]) +
    theme_std(plotDefsObj) +
    theme(
      ## axis.text.y = element_blank(),
    ) +
    guides(
      group = expObjFactors[["factorList"]][["group"]]
    )
    ## scale_fill_identity(
    ##   defs["names.heatmap.geneCounts.legend"],
    ##   labels = factor(cats, levels = rev(cats)),
    ##   breaks = factor(breaks),
    ##   guide = "legend"
    ## )

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

getAnalyses <- function() {
  conf <- getConf()
  analyses <- list()
  if (conf$analyses$go_tf == TRUE) {
    analyses <- append(analyses, list(c(
                                   go_regulation_of_dna_transcription_0006355 = "transcription_reg",
                                   ptrpf2 = "ptrpf2",
                                   pf2hyp1 = "ptrpf2_like")))
  }

  if (conf$analyses$go_cwde == TRUE) {
    analyses <- append(analyses, list(c(cwde = "cwde")))
  }

  if (conf$analyses$go_transport == TRUE) {
    analyses <- append(analyses, list(c(go_transmembrane_transport_0055085 = "transmembrane_transport")))
  }

  if (conf$analyses$effectors == TRUE) {
    analyses <- append(analyses, list(c(
                                   effectors = "effectors",
                                   toxa = "toxa")))
  }

  if (conf$analyses$go_kinases) {
    analyses <- append(analyses, list(c(
                                    go_kinase_serine_threonine_0004674 = "ser_thr_kinases"
                                  )))
  }

  if (conf$analyses$go_intracellular_signalling) {
    analyses <- append(analyses, list(c(
                                    go_intracellular_signalling_0035556 = "intracellular_signalling"
                                  )))
  }

  if (conf$analyses$go_chromatin) {
    analyses <- append(analyses, list(c(
                                   go_chromatin_remodeling_0006338 = "chromatin_remodeling",
                                   chromatin_remodeling_factors = "chromatin_remodeling_factors",
                                   jmjc_histone_demethylases = "jmjc_demethylases",
                                   histone_acetyltransferases = "acetyltransferases",
                                   histone_methyltransferases = "methyltransferases"
                                  )))
  }

  if (conf$analyses$go_methylation) {
    analyses <- append(analyses, list(c(
                                   go_histone_methyltransferase_0042054 = "methyltransferase",
                                   dcm_dna_methyltransferases = "dna-cytosine_methyltransferase"
                                  )))
  }

  ## if (conf$analyses$go_rnai) {
  ##   analyses <- append(analyses, list(c(
  ##                                 )))
  ## }

  ## if (conf$analyses$go_transposons) {
  ##   analyses <- append(analyses, list(c(
  ##                                 )))
  ## }

  return(analyses)
}

goHeatmaps <- function(counts, zscore, clustering = TRUE) {
  if (zscore == TRUE) {
    counts <- (scale(counts))
  }

  analyses <- getAnalyses()
  plotList <- list()
  for (a in 1:length(analyses)) {
    d1 <- new.env()

    ## TODO: Remove duplicated genes from any single category and move
    ## genes that appear in more than one category into a unique one
    d1$cats <- analyses[[a]]

    d1$genes <-
      (function(x, y) {
        sig <<- x[names(x) == names(y)[1]][x[names(x) == names(y)[1]] %in% sdeg$names]
        if (length(x[names(x) == names(y)[1]]) > 100) {
          ret <- sig
        } else {
          ret <- x[names(x) == names(y)[1]]
        }
        return(ret)
      })(geneList$genes, d1$cats)
    d1$counts <- counts[d1$genes, ]

    if (length(d1$cats) > 1) {
      for (i in 2:length(d1$cats)) {
        d1$tmp <- counts[geneList$genes[names(geneList$genes) == names(d1$cats)[i]], ]

        if (d1$tmp %>% is.vector()) {
          d1$tmp <- t(data.frame(d1$tmp))
          rownames(d1$tmp) <- geneList$genes[names(geneList$genes) == names(analyses[[a]])[i]]
        }

        d1$genes <- append(d1$genes, geneList$genes[names(geneList$genes) == names(analyses[[a]])[i]])
        d1$counts <- rbind(d1$counts, d1$tmp)
      }

      for (i in 1:length(d1$cats)) {
        d1$genes[names(d1$genes) == names(d1$cats)[i]] <-
          paste(
            d1$cats[i],
            (d1$genes[names(d1$genes) == names(d1$cats)[i]]),
            sep = ":"
          )
        ## (function(x, y) rownames(x)[rownames(x) == (y[names(y) == names(d1$cats)[i]])]) (d1$counts, d1$genes), sep = ":")
      }
    }

    ## paste(
    ##   d1$cats[i],
    ##   (d1$genes[names(d1$genes) == names(d1$cats)[i]]),
    ##   sep = ":"
    ## )

    if (d1$counts %>% is.vector()) {
      d1$counts <- t(data.frame(d1$counts))
    }
    rownames(d1$counts) <- d1$genes

    d1$remove <- c()
    for (i in 1:nrow(d1$counts))
      for (j in 1:ncol(d1$counts))
        if (is.na(d1$counts[i, j])) {d1$remove <- append(d1$remove, i); break}
    if (length(d1$remove) > 0) d1$counts <- d1$counts[-d1$remove, ]

    d1$remove <- grep(".*[.]1.*", rownames(d1$counts))
    if (length(d1$remove) > 0) d1$counts <- d1$counts[-d1$remove, ]

    for (i in 1:nrow(d1$counts)) {
      for (g in sig) {
        d1$tmp <- grep(qq(".*@{g}.*"), rownames(d1$counts)[i])
        if (length(d1$tmp) > 0) {rownames(d1$counts)[i] <- paste0(rownames(d1$counts)[i], "*"); break}
      }
    }

    show_rows <- TRUE
    if (length(d1$genes) > 100) {show_rows <- FALSE}

    colours <- list(
      cultivar = hcl.colors(length(expFactors$factorList$cultivar), "Set 3"),
      time = hcl.colors(length(expFactors$factorList$time), "Dark 3"))

    names(colours$cultivar) <- expFactors$factorList$cultivar
    names(colours$time) <- expFactors$factorList$time

    d1$plots <- list()
    for (t in 1:length(expFactors$factorList$time)) {
      d1$tcounts <-
        d1$counts[, paste0(expFactors$factorList$cultivar, toTitleCase(expFactors$factorList$time[[t]]))]
      if (d1$tcounts %>% is.vector) d1$tcounts <- t(data.frame(d1$tcounts))
      rownames(d1$tcounts) <- rownames(d1$counts)

      ## if (clustering == TRUE) {clusters <- c(TRUE, (if (t == 1) {TRUE} else {FALSE}))}
      clusters <- c(TRUE, TRUE)
      d1$plots[[t]] <- pheatmap(
        as.matrix(d1$tcounts),
        main = (paste(analyses[[a]][[1]], expFactors$factorList$time[[t]])),
        color = rev(hcl.colors(3, "RdYlGn")),
        show_rownames = show_rows,
        cluster_cols = clusters[1], cluster_rows = clusters[2],
        annotation_col = data.frame(
          ## groups = data.frame(expFactors$factorList$group),
          cultivar =
            rep(expFactors$factorList$cultivar, each = expData$nfactors - 1)
          ## time = rep(expFactors$factorList$time, length(expFactors$factorList$group) / length(expFactors$factorList$time))
        ),
        annotation_colors = colours
      )
    }

    d1$plots[[t + 1]] <- pheatmap(
      as.matrix(d1$counts),
      main = (paste(analyses[[a]][[1]], "all time points", sep = ", ")),
      color = rev(hcl.colors(3, "RdYlGn")),
      show_rownames = show_rows,
      cluster_cols = TRUE, cluster_rows = TRUE,
      annotation_col = data.frame(
        ## groups = data.frame(expFactors$factorList$group),
        cultivar = rep(expFactors$factorList$cultivar, each = expData$nfactors),
        time = rep(expFactors$factorList$time, length(expFactors$factorList$group) / length(expFactors$factorList$time))
      ),
      annotation_colors = colours
    )

    plotList[[a]] <- d1$plots
  }

  return(plotList)
}

saveHeatmaps <- function(list = heatmapList) {
  analyses <- getAnalyses()
  for (a in 1:length(analyses)) {
    for (t in 1:(length(expFactors$factorList$time) + 1)) {
      if (t == length(expFactors$factorList$time) + 1) {
        tstr = "all"
      } else {
        tstr = expFactors$factorList$time[t]
      }

      f = formatPath(conf$data$plots, species, "heatmaps", tstr, root = FALSE)
      p <- list[[a]][[t]]
      h = round(200 + 20 * round(dim(p)[1], 10), 100)
      w = round(100 + 10 * round(10 * dim(p)[2], 10), 100)

      makePath(f)
      ## png(qq("@{f}/@{analyses[[a]][1]}.png"), height = h, width = w)
      pdf(file = qq("@{f}/@{analyses[[a]][1]}.pdf"), paper = "a4")
      plot(p)
      dev.off()
    }
  }
}
