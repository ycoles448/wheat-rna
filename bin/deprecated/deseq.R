getDifferentialGenes <- function(resultsList, dds, expObj, expObjFactors, alpha = SIG, l2fc = L2FC) {
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

plotDifferentialGenes <- function(sdegList, expObjFactors, plotDefsObj, offsetR = 0, offsetC = 0) {
  k <- length(expObjFactors[["factorList"]][["group"]]) - 1

  nfactors <- length(factors)
  ngroups <- (k + 1) / nfactors

  plotList <- list()
  # Multiple Venn diagrams are too messy...
  # Combining different groups into a Venn diagram!
  # Construct pairs by row
  for (row in nfactors * 1:(ngroups) - 1 + offsetR) {
    rowname <- names(sdegList[[length(sdegList)]])[row]
    colname <- names(sdegList)[1]
    for (col in nfactors * 1:(ngroups - 1) - 1 + offsetC) {
      colname <- names(sdegList)[col]

      if (colname == rowname) {
        col <- length(sdegList) - 1
        colname <- names(sdegList)[col]
      }

      printf("row: %s\t rowname: %s\tcol: %s\tcolname: %s\n", row, rowname, col, colname)
      ## sdegList[[row]][[col]]
    }
    printf("\n")
  }


  ## # Iterate over number of factors
  ## for (i in nfactors:1) {
  ##   # Iterate over number of groups
  ##   ## printf("Iter: %s\n", i)
  ##   i2 <- (nfactors:1)[i]

  ##   plotList[[i2]] <- list()
  ##   for (j in nfactors * (1:(ngroups))) {
  ##     j2 <- j / nfactors
  ##     col <- nfactors - i + 1
  ##     n <- 1

  ##     plotList[[i2]][[j2]] <- list()
  ##     for (row in nfactors * (1:(ngroups - 1)) + (nfactors:1)[i] - 1) {
  ##       printf("Row: %s\t Col: %s\t Group: %s\n", row, col, j)

  ##       contrast <- sdegList[[row]][[col]][["contrast"]]
  ##       printf("%s\n\n", contrast)
  ##       plotList[[i2]][[j2]][["up"]][[contrast]] <- sdegList[[row]][[col]][["up"]]
  ##       plotList[[i2]][[j2]][["dn"]][[contrast]] <- sdegList[[row]][[col]][["dn"]]

  ##       if (col < j - 1) {
  ##         col <- col + nfactors
  ##       }

  ##       n <- n + 1
  ##     }
  ##   }
  ## }

  return(plotList)

  ## # Create an empty plot
  ## pEmpty <- ggVennDiagram(
  ##   list(0, 0),
  ##   label_alpha = 0
  ## ) +
  ##   theme_std(plotDefsObj) +
  ##   theme_empty(plotDefsObj) +
  ##   guides_none()

  ## for (r in 1:k) {
  ##   # Fill with data
  ##   for (c in 1:n) {
  ##     ## printf("r: %s\tc: %s\n", r, c)
  ##     p <- ggVennDiagram(sdegList[[r]][[c]], label_size = 2) +
  ##       theme_std(plotDefsObj) +
  ##       theme_venn(plotDefsObj) +
  ##       guides_none()

  ##     plotList[[m]] <- ggplotGrob(p)
  ##     m <- m + 1
  ##   }

  ##   # Fill with empty plots
  ##   n <- n + 1
  ##   for (c in n:k) {
  ##     plotList[[m]] <- ggplotGrob(pEmpty)
  ##     m <- m + 1
  ##   }
  ## }

  ## for (r in 1:k) {
  ##   pr <- plotList[[k * (r - 1) + 1]]
  ##   for (c in 2:k) {
  ##     i <- k * (r - 1) + c
  ##     ## printf("Index: %s\n", i)
  ##     pr <- rbind(plotList[[i]], pr)
  ##   }

  ##   if (!exists("plotArr")) {
  ##     plotArr <- pr
  ##   } else {
  ##     plotArr <- cbind(plotArr, pr)
  ##   }
  ## }

  ## return(plotArr)
}
