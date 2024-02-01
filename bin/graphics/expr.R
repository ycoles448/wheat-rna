plotExprProfilesGo <- function(clustPlots, goMap, globalColours = FALSE) {
  goList <- list()
  goTermList <- list()
  for (j in 1:length(modules)) {
    goList[[j]] <- list()
    clustPlots[[j]] <- list()
  }

  # Get GO-term enrichment
  for (i in 1:length(ont)) {
    for (j in 1:length(modules)) {
      genes <- modules[[j]][which(modules[[j]] %in% rownames(countsRL))]
      genesGo <- factor(rep(1, length = length(genes)), levels = c(0, 1))
      names(genesGo) <- genes

      go <- new(
        "topGOdata",
        ontology = ont[i],
        allGenes = genesGo,
        annotationFun = annFUN.gene2GO,
        gene2GO = goMap,
      )

      goTest <- runTest(go, algorithm = "classic", statistic = "ks")
      goTestTable <- GenTable(go, classicKS = goTest, orderBy = "classicKS")

      goList[[j]][[i]] <- goTestTable
    }
  }

  # Determine universe of GO-term annotations
  for (i in 1:length(ont)) {
    goTermList[[i]] <- vector("character")
    for (j in 1:length(modules)) {
      goTermList[[i]] <- append(goTermList[[i]], goList[[j]][[i]][, "Term"])
    }
    goTermList[[i]] <- unique(goTermList[[i]])
  }

  # Universal colour mapping of terms
  goTermColours <- list()
  for (i in 1:length(goTermList)) {
    # Too many colours for viridis colour sets
    ## switch(i,
    ##   {
    ##     colFun <- magma
    ##   },
    ##   {
    ##     colFun <- mako
    ##   },
    ##   {
    ##     colFun <- rocket
    ##   },
    ## )
    colFun <- rainbow

    goTermColours[[i]] <- colFun(length(goTermList[[i]]))
    names(goTermColours[[i]]) <- goTermList[[i]]
  }

  # Distribution bars
  for (i in 1:length(ont)) {
    for (j in 1:length(modules)) {
      d1 <- data.frame(
        "term" = goList[[j]][[i]][, "Term"],
        "count" = goList[[j]][[i]][, "Annotated"],
        "group" = rep(1, nrow(goList[[j]][[i]])),
        "col" = rep(NA, nrow(goList[[j]][[i]]))
      )

      # Fix bar order using random
      rand <- sort(runif(length(1e3 * d1[, "count"]), 0, 1))
      for (k in 1:nrow(d1)) {
        d1[k, "count"] <- d1[k, "count"] + rand[k]
      }
      d1 <- d1[order(d1[, "count"]), ]
      d1[, "count"] <- factor(d1[, "count"], levels = d1[, "count"])

      for (k in 1:nrow(d1)) {
        d1[k, "col"] <- goTermColours[[i]][which(d1[k, "term"] == names(goTermColours[[i]]))]
      }

      p <- ggplot(d1, aes(x = count, y = group, fill = col)) +
        geom_bar(stat = "identity", position = "stack", orientation = "y") +
        theme_std_highres(plotDefs) +
        theme(
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
        ) +
        guides(fill = guide_legend(nrow = 4, title = qq("Ontology: @{ont[[i]]}"))) +
        xlab("Distribution of GO terms")

      if (globalColours == TRUE) {
        p <- p + scale_fill_identity(
          breaks = rev(d1[, "col"]),
          labels = rev(d1[, "term"]),
          guide = "legend"
        )
      } else {
        p <- p + scale_fill_brewer(
          name = qq("Term: @{ont[[i]]}"),
          palette = "Set3",
          breaks = rev(d1[, "col"]),
          labels = rev(d1[, "term"]),
        )
      }

      pLeg <- cowplot::get_legend(p + theme(legend.position = NULL))

      clustPlots[[j]][[i]] <- list(p, pLeg)
    }
  }

  return(clustPlots)
}

plotExprProfilesGoRestricted <- function(clustPlots, goMap, globalColours = FALSE, db) {
  goTermList <- list()
  for (i in 1:length(ont)) {
    if (ont[i] == "BP") {
      fileGoNames <- formatPath(conf[["data"]][["data"]], conf[["data"]][["go"]], species, qq("biological_process-@{db}.txt"))
    } else if (ont[i] == "MF") {
      fileGoNames <- formatPath(conf[["data"]][["data"]], conf[["data"]][["go"]], species, qq("molecular_function-@{db}.txt"))
    } else if (ont[i] == "CC") {
      fileGoNames <- formatPath(conf[["data"]][["data"]], conf[["data"]][["go"]], species, qq("cellular_component-@{db}.txt"))
    }

    print(qq("Database file: @{fileGoNames}"))

    goTermList[[i]] <- vector("character")
    for (j in 1:length(modules)) {
      goNames <- read.csv(fileGoNames, header = FALSE, sep = "\t")
      d1 <- goNames[which(goNames[, 1] %in% modules[[j]]), 2]
      goTermList[[i]] <- append(goTermList[[i]], unlist(sapply(d1, function(x) {
        strsplit(x, ";")
      })))
    }

    goTermList[[i]] <- unique(goTermList[[i]])
  }

  goTermColours <- list()
  for (i in 1:length(goTermList)) {
    switch(i,
      {
        colFun <- viridis
      },
      {
        colFun <- plasma
      },
      {
        colFun <- turbo
      },
    )
    goTermColours[[i]] <- colFun(length(goTermList[[i]]))
    names(goTermColours[[i]]) <- goTermList[[i]]
  }

  for (i in 1:length(modules)) {
    clustPlots[[i]] <- list()
    for (j in 1:length(ont)) {
      if (ont[j] == "BP") {
        fileGoNames <- formatPath(conf[["data"]][["data"]], conf[["data"]][["go"]], species, qq("biological_process-@{db}.txt"))
      } else if (ont[j] == "MF") {
        fileGoNames <- formatPath(conf[["data"]][["data"]], conf[["data"]][["go"]], species, qq("molecular_function-@{db}.txt"))
      } else if (ont[j] == "CC") {
        fileGoNames <- formatPath(conf[["data"]][["data"]], conf[["data"]][["go"]], species, qq("cellular_component-@{db}.txt"))
      }

      goNames <- read.csv(fileGoNames, header = FALSE, sep = "\t")
      goColours <- goTermColours[[j]]

      d1 <- goNames[which(goNames[, 1] %in% modules[[i]]), ]
      terms <- unlist(sapply(d1[, 2], function(x) {
        unlist(strsplit(x, ";"))
      }))

      ## goCount <- vector("integer", length = length(unique(terms)))
      ## goCol <- vector("integer", length = length(unique(terms)))
      goCount <- c()
      goCol <- c()
      for (k in 1:length(unique(terms))) {
        goCount[k] <- sum(unique(terms)[k] == terms)
        goCol[k] <- goColours[unique(terms)[k] == names(goColours)]
      }

      # Fix bar order using random
      goCount <- 1e3 * goCount
      rand <- sort(runif(length(goCount), 0, 1))
      for (k in 1:length(goCount)) {
        goCount[k] <- goCount[k] + rand[k]
      }

      d2 <- data.frame(
        group = rep(1, length(unique(terms))),
        cat = as.character(unique(terms)),
        col = as.character(goCol),
        count = factor(
          as.numeric(goCount),
          levels = sort(unique(goCount, decreasing = FALSE))
        )
        ## count = as.integer(goCount)
      )
      d2 <- d2[order(d2[, "count"]), ]
      ## d2 <- d2 %>% arrange(count)

      p <- ggplot(d2, aes(x = count, y = group, fill = col)) +
        geom_bar(stat = "identity", position = "stack", orientation = "y") +
        ## scale_fill_brewer(palette = "Set3") +
        ## scale_fill_viridis_d(option = "H") +
        ## scale_fill_identity(
        ##   breaks = rev(d2[, "col"]),
        ##   labels = rev(d2[, "cat"]),
        ##   guide = "legend"
        ## ) +
        theme_std_highres(plotDefs) +
        theme(
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
        ) +
        guides(
          fill = guide_legend(nrow = 4, title = paste(ont[j], "(Aspergillus)"))
        ) +
        xlab("Distribution of GO terms")

      if (globalColours == TRUE) {
        p <- p + scale_fill_identity(
          breaks = rev(d2[, "col"]),
          labels = rev(d2[, "cat"]),
          guide = "legend"
        )
      } else {
        ## p <- p + scale_fill_brewer(
        ##   name = qq("Term: @{ont[[i]]}"),
        ##   palette = "Set3",
        ##   breaks = rev(d2[, "col"]),
        ##   labels = rev(d2[, "cat"]),
        ## )
        p <- p + scale_fill_viridis_d(
          option = "H",
          name = qq("Term: @{ont[[i]]}"),
          breaks = rev(d2[, "col"]),
          labels = rev(d2[, "cat"]),
          )
      }

      pLeg <- cowplot::get_legend(p + theme(legend.position = NULL))

      clustPlots[[i]][[j]] <- list(p, pLeg)
    }
  }

  return(clustPlots)
}
