loadData <- function(mat, meta, contrast, factors, nreps, species) {
  expData <- expObj(
    mat = mat,
    meta = meta,
    contrast = as.character(contrast),
    factors = as.vector(factors),
    nfactors = as.integer(length(factors)),
    nreps = as.integer(nreps),
    species = as.character(species)
  )

  return(expData)
}

fixData <- function(expObj) {
  colnames(expObj[["mat"]]) <- rownames(expObj[["meta"]])
}

loadExperimentData <- function(expObj, geneList, threshold) {
  toRm <- c()
  countMat <- expObj[["mat"]]
  for (i in 1:nrow(countMat)) {
    # Remove genes, if sum of counts across all samples does not exceed threshold
    ## if (sum(countMat[i, ]) < threshold && !(rownames(countMat)[i] %in% geneList[["genes"]])) {
    ##   toRm <- append(toRm, i)
    ## }

    # Remove genes, if less than 50% of samples exceed the threshold
    c <- 0
    if (sum(countMat[i, ] > threshold, na.rm = TRUE) > (ncol(countMat) / 2) && !(rownames(countMat)[i] %in% geneList[["genes"]])) {
      toRm <- append(toRm, i)
    }
  }

  countMat <- countMat[-toRm, ]

  x <- DESeqDataSetFromMatrix(
    countData = countMat,
    colData = expObj[["meta"]],
    design = formula(paste0("~", expObj[["contrast"]]))
  )
  dds <- DESeq(x, parallel = TRUE)

  return(dds)
}

loadGeneList <- function(path) {
  # Return a named vector with k = category, v = gene
  sep <- "/" # TODO: Change file separator, based on OS
  absPath <- paste(getwd(), path, sep = sep)
  files <- list.files(absPath)

  # Also TODO: use formatPath()

  genes <- c()
  for (f in files) {
    name <- file_path_sans_ext(f)
    buffer <- fread(paste(absPath, f, sep = sep), header = FALSE)
    x <- 0
    for (l in 1:length(buffer[[1]])) {
      x <- x + 1
      genes[x] <- buffer[[1]][l]
      names(genes)[x] <- name
    }
  }

  return(geneListObj(genes = genes))
}

addSelectGenes <- function(geneListObj, genesSelect) {
  len <- length(geneListObj[["genes"]])
  for (i in 1:length(genesSelect)) {
    geneListObj[["genes"]][len + i] <- genesSelect[i]
    names(geneListObj[["genes"]])[len + i] <- names(genesSelect)[i]
  }
}

getTotalCounts <- function(counts, n) {
  # Get the average of total counts per gene
  totalCounts <- c()
  for (i in n * (1:(ncol(counts) / n))) {
    totalCounts[i / n] <- sum(counts[, (i - n + 1):i]) / n
  }

  return(totalCounts)
}

stripNumbers <- function(x) {
  c <- 0
  lastChar <- substr(x, nchar(x) - c, nchar(x) - c)
  num <- lastChar

  while (!is.na(as.numeric(lastChar))) {
    c <- c + 1
    lastChar <- substr(x, nchar(x) - c, nchar(x) - c)
    num <- paste0(lastChar, num)
  }

  num <- substr(num, 2, nchar(num)) # Remove trailing letter
  return(c(substr(x, 1, nchar(x) - c), num))
}

iterDifferentialComplex <- function(x, argList, expObjFactors, offsetR, offsetC) {
  k <- length(expObjFactors[["factorList"]][["group"]]) - 1
  factors <- expObjFactors[["factorList"]][["factors"]]
  nfactors <- length(factors)

  for (i in nfactors * 1:((k + 1) / 2) - 1 + offsetR) {
    c <- nfactors + offsetC
    for (k in nfactors * 1:((k + 1) / 2) - 1 + offsetR) {
      r <- j
      if (r == nfactors && c == nfactors) {
        next
      }

      while (r < i) {
        r <- r + nfactors
      }

      do.call(x, argList)

      if (c < i) {
        c <- c + nfactors
      }
    }
  }
}

getModules <- function(clustRuns) {
  counter <- c()
  for (i in 1:length(clustRuns)) {
    clustRuns[[i]] <- read.table(
      formatPath(dirClust, names(clustRuns)[i], "Clusters_Objects.tsv"),
      header = TRUE, sep = "\t"
    )
    counter <- append(counter, dim(clustRuns[[i]])[2])
  }

  # Only use Clust run with largest number of distinct modules
  modules <- clustRuns[[which(counter == max(counter))]]

  return(modules)
}
