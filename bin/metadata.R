getFactors <- function(factorList, factors) {
  f <- list()
  for (i in 1:length(factors)) {
    name <- factors[i]
    nameOrdered <- paste0(name, "Ordered")

    f[[name]] <- factorList[[name]]

    if (!is.null(factorList[[nameOrdered]])) {
      f[[nameOrdered]] <- factorList[[nameOrdered]]
    } else {
      f[[nameOrdered]] <- factorList[[name]]
    }
  }

  expFactors <- expObjFactors(
    factorList = f
  )

  return(expFactors)
}

setFactorSamples <- function(expObj, expObjFactors, group) {
  g <- expObj[["meta"]][[expObj[["contrast"]]]]
  s <- rownames(expObj[["meta"]])

  for (i in 1:length(g)) {
    expObjFactors[["factorList"]][["sample"]][i] <- paste0(g[i], s[i])
  }

  # Set ordered factors for sample and group
  # Get factor names from contrast
  n <- strsplit(tolower(
    gsub(
      "(?!^)([[:upper:]])",
      " \\1",
      substr(expObj[["contrast"]], nchar("group") + 1, nchar(expObj[["contrast"]])),
      perl = TRUE
    )
  ), " ")[[1]]

  # Get factors strictly from names
  f <- list()
  for (i in 1:length(n)) {
    f[[n[[i]]]] <- expObjFactors[["factorList"]][[paste0(n[i], "Ordered")]]
  }

  # Get linear combinations of factors
  groupOrdered <- c()
  for (i in 1:length(f)) {
    if (i == length(f)) {
      break
    }

    for (j in f[[i]]) {
      for (k in f[[i + 1]]) {
        groupOrdered <- append(groupOrdered, paste0(toTitleCase(j), toTitleCase(k)))
      }
    }
  }

  # Map the unordered sample to their respective ordered group
  x <- expObjFactors[["factorList"]][["sample"]] # Unordered sample name
  y <- groupOrdered # Ordered group name
  z <- vector(length = length(y))
  for (i in 1:length(x)) {
    group <- stripNumbers(x[i])[1]
    num <- as.integer(stripNumbers(x[i])[2])

    for (j in 1:length(y)) {
      if (group == y[j]) {
        ## printf(paste0("LHS (x): ", group, "\t\tRHS (y): ", y[j], "\t\tx index: ", i, "\t\ty index: ", j, "\t\tNum: ", num, "\n"))
        z[j] <- i
        break
      }
    }
  }

  # Assign each sample in order
  k <- 0
  for (i in z) {
    for (j in expObj[["nreps"]]:1) {
      k <- k + 1
      expObjFactors[["factorList"]][["sampleOrdered"]][k] <- x[(i - j + 1)]
    }
  }

  for (i in 1:length(n)) {

  }

  f <- list()
  for (i in 1:length(n)) {
    f[[n[[i]]]] <- expObjFactors[["factorList"]][[n[i]]]
  }

  # Get linear combinations of factors, unordered
  group <- c()
  for (i in 1:length(f)) {
    if (i == length(f)) {
      break
    }

    for (j in f[[i]]) {
      for (k in f[[i + 1]]) {
        group <- append(group, paste0(toTitleCase(j), toTitleCase(k)))
      }
    }
  }

  expObjFactors[["factorList"]][["group"]] <- group
  expObjFactors[["factorList"]][["groupOrdered"]] <- groupOrdered
}

setFactor <- function(expObj, expObjFactors, list) {
  for (i in 1:length(list)) {
    expObjFactors[["factorList"]][[names(list)[i]]] <- list[[i]]
  }
}

exportToClust <- function(filename, ref, group, meta) {
  # Remove old file
  unlink(filename)

  d <- data.frame(
    group = meta[, group]
  )
  rownames(d) <- meta[, "sample"]

  # Write Clust data
  f <- file(filename)
  text <- ""
  for (g in unique(d[, "group"])) {
    reps <- paste(rownames(d)[d[, "group"] == g], collapse = ",")
    text <- paste0(text, qq("@{ref}\t@{g}\t@{reps}\n"))
  }
  writeLines(text, f)
  close(f)
}

exportTraits <- function() {
  # Generic function to export traits based off of factors
}

exportTraitsFungi <- function(filename, meta, factors, SEP = "\t") {
  unlink(filename)
  f <- file(filename)
  text <- ""

  # Get unique values for each factor
  factorValues <- list()
  for (i in factors) {
    factorValues[[i]] <- unique(meta[, i])
  }

  # Create list of traits from all factors
  traits <- c()
  for (i in 1:length(factorValues)) {
    for (t in factorValues[[i]]) {
      traits <- append(traits, paste0(names(factorValues)[i], ":", t))
    }
  }

  # Header
  text <- paste0(text, "sample", SEP)
  for (i in 1:length(traits)) {
    text <- paste0(text, traits[i], SEP)
  }
  text <- paste0(text, "\n")

  for (i in 1:nrow(meta)) {
    text <- paste0(text, rownames(meta)[i], SEP)
    for (j in 1:length(traits)) {
      # Check if matrix value is one of the cultivars
      t <- traits[j]
      col <- strsplit(t, ":")[[1]][1]
      val <- strsplit(t, ":")[[1]][2]
      if (meta[i, col] == val) {
        check <- 1
      } else {
        check <- -1
      }
      text <- paste0(text, check, SEP)
    }
    text <- substr(text, 1, nchar(text) - 1) # Remove trailing separator
    text <- paste0(text, "\n")
  }

  writeLines(text, f)
  close(f)
}
