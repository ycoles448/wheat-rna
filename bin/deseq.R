getResultsList <- function(dds, expObj, expObjFactors, alpha = SIG) {
  k <- length(expObjFactors[["factorList"]][["group"]]) - 1
  n <- 1

  resultsList <- list()
  for (r in 1 + 1:k) {
    gr <- expObjFactors[["factorList"]][["group"]][r]
    resultsList[[gr]] <- list()

    for (c in 1:n) {
      ## printf("r:%s\tc:%s\n", r, c)
      ## printf("gr:%s\tgc:%s\n\n", gr, gc)
      gc <- expObjFactors[["factorList"]][["group"]][c]

      resultsList[[gr]][[gc]] <- results(
        dds,
        alpha = alpha,
        contrast = c(expObj[["contrast"]], gr, gc),
        parallel = TRUE
      )
    }

    n <- n + 1
  }

  return(resultsList)
}

getDifferentialGenes <- function(resultsList, expObjFactors, alpha = SIG, l2fc = L2FC) {
  k <- length(expObjFactors[["factorList"]][["group"]])
  n <- 1

  sdegList <- list()
  for (r in 1:k) {
    gr <- expObjFactors[["factorList"]][["group"]][r]
    sdegList[[gr]] <- list()

    for (c in 1:n) {
      ## printf("r: %s\tc: %s\n", r, c)
      gc <- expObjFactors[["factorList"]][["group"]][c]

      if (gr == gc) {
        resultsList[[gr]][[gc]] <- NA
        next
      }

      mat <- resultsList[[gr]][[gc]]
      mat <- subset(mat, !is.na(padj))
      contrast <- paste(strsplit(attr(mat, "elementMetadata")[2, 2], " ")[[1]][6:8], collapse = " ")
      up <- mat[mat[, "padj"] < SIG & mat[, "log2FoldChange"] > l2fc, ]
      dn <- mat[mat[, "padj"] < SIG & mat[, "log2FoldChange"] < l2fc, ]
      sdegList[[gr]][[gc]] <- list(contrast = contrast, up = rownames(up), dn = rownames(dn))
    }

    n <- n + 1
  }

  return(sdegList)
}

plotDifferentialGenes <- function(sdegList, expObjFactors, plotDefs, offsetR = 0, offsetC = 0, debug = FALSE) {
  # Venn diagram
  d <- plotDefs[["defs"]]
  k <- length(expObjFactors[["factorList"]][["group"]]) - 1
  factors <- expObjFactors[["factorList"]][["factors"]]
  nfactors <- length(factors)

  # Get colour limits
  if (debug == FALSE) {
    upCount <- c()
    dnCount <- c()
    for (i in nfactors * 1:((k + 1) / 2) - 1 + offsetR) {
      group <- names(sdegList)[i + 1]
      c <- nfactors + offsetC

      genesUp <- list()
      genesDn <- list()
      for (j in nfactors * 1:((k + 1) / 2) + offsetR) {
        r <- j
        if (r == nfactors && c == nfactors) {
          next
        }

        while (r < i) {
          r <- r + nfactors
        }

        upCount <- append(upCount, length(sdegList[[r]][[c]][["up"]]))
        dnCount <- append(dnCount, length(sdegList[[r]][[c]][["dn"]]))

        if (c < i) {
          c <- c + nfactors
        }
      }
    }

    minUp <- min(upCount)
    maxUp <- max(upCount)

    minDn <- min(dnCount)
    maxDn <- max(dnCount)

    printf("minUp: %s\tmaxUp: %s\tminDn: %s\tmaxDn: %s\n", minUp, maxUp, minDn, maxDn)
    ## printf("\n")
  }

  plotList <- list()
  for (i in nfactors * 1:((k + 1) / 2) - 1 + offsetR) {
    group <- names(sdegList)[i + 1]
    c <- nfactors + offsetC

    genesUp <- list()
    genesDn <- list()
    catNames <- c()
    for (j in nfactors * 1:((k + 1) / 2) + offsetR) {
      r <- j
      if (r == nfactors && c == nfactors) {
        next
      }

      while (r < i) {
        r <- r + nfactors
      }

      # Do stuff here
      ## printf("r: %s\tc: %s\n", r, c)
      rname <- names(sdegList)[r]
      cname <- names(sdegList[[r]])[c]
      contrast <- sdegList[[r]][[c]][["contrast"]]
      ## printf("i: %d\tj: %s\tr: %s\tc: %s\trname: %s\tcname: %s\tcontrast: %s\n", i, j, r, c, rname, cname, contrast)

      genesUp[[contrast]] <- sdegList[[r]][[c]][["up"]]
      genesDn[[contrast]] <- sdegList[[r]][[c]][["dn"]]

      if (c == i) {
        ## printf("Col 2 selected\n");
        catNames <- append(catNames, strsplit(contrast, " ")[[1]][1])
      } else {
        ## printf("Col 1 selected\n");
        catNames <- append(catNames, strsplit(contrast, " ")[[1]][3])
      }

      # End of stuff
      if (c < i) {
        c <- c + nfactors
      }
    }
    ## print(str(sdegList[[r]]))

    if (debug == FALSE) {
      # Set category names
      ## m <- as.numeric(d["names.venn.margin"])
      catNames <- sapply(names(genesUp), function(x) {
        ## x <- gsub(" vs ", " vs\n", x)
        names <- strsplit(x, " vs ")[[1]]
        x <- names[which(names != strsplit(group, " ")[[1]][1])]
        return(x)
      })
      ## print(catNames)
      ggVennDiagramStd <- function(data, catNames) {
        return(ggVennDiagram(
          data,
          catNames,
          label = "count",
          set_size = 3,
          label_alpha = 0.5,
          label_size = 3,
          label_geom = "label",
        ) +
          ## theme_minimal() +
          theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
          scale_x_continuous(expand = expansion(0.1)))
      }

      pUp <- ggVennDiagramStd(genesUp, catNames) +
        ggtitle(toTitleCase(strsplit(group, " ")[[1]][1])) +
        scale_fill_continuous(
          limits = c(minUp, maxUp)
        ) +
        guides_none()
      ## theme_venn(plotDefs)
      pDn <- ggVennDiagramStd(genesDn, catNames) +
        ggtitle(toTitleCase(strsplit(group, " ")[[1]][1])) +
        scale_fill_continuous(
          limits = c(minDn, maxDn)
        ) +
        guides_none()
      ## theme_venn(plotDefs)

      plotList[[group]] <- list()
      plotList[[group]][["up"]] <- pUp
      plotList[[group]][["dn"]] <- pDn

      ## Extract the legend from a separate plot
      pLegUp <- get_legend(
        ggVennDiagram(genesUp) +
          scale_fill_continuous(limits = c(minUp, maxUp))
      )
      pLegDn <- get_legend(
        ggVennDiagram(genesDn) +
          scale_fill_continuous(limits = c(minDn, maxDn))
      )

      legends <- list(up = pLegUp, dn = pLegDn)
    }
  }

  plotList <- list(plotList, legends)
  return(plotList)
}
