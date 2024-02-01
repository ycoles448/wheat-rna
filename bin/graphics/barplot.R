plotTotalCountsBars <- function(counts, expObjFactors, plotDefsObj) {
  d <- reshape2::melt(counts)
  defs <- plotDefsObj[["defs"]]
  colnames(d) <- "counts"

  d[, "group"] <- factor(
    expObjFactors[["factorList"]][["group"]],
    levels = expObjFactors[["factorList"]][["groupOrdered"]]
  )

  p <- ggplot(d, aes(x = group, y = counts, fill = group)) +
    geom_col(color = "black") +
    ggtitle(defs["names.barplot.totalCounts.title"]) +
    xlab(defs["names.xlab.cultivar"]) +
    ylab(defs["names.ylab.totalCounts"]) +
    theme_std(plotDefsObj) +
    guides_none()

  return(p)
}

plotDiffExprTime <- function(sdegGroups, sdeg) {
  getDF1 <- function(x, f1, f2) {
    d[paste0(sprintf(fmtstr, f1), sprintf(fmtstr, f2)), "dn"] <<- -count(x[, "log2FoldChange"] < -SIG)
    d[paste0(sprintf(fmtstr, f1), sprintf(fmtstr, f2)), "up"] <<- +count(x[, "log2FoldChange"] > +SIG)
    d[paste0(sprintf(fmtstr, f1), sprintf(fmtstr, f2)), "dnannot"] <<-
      -count(rownames(x[x[, "log2FoldChange"] < -SIG, ]) %in% goMap[, 1])
    d[paste0(sprintf(fmtstr, f1), sprintf(fmtstr, f2)), "upannot"] <<-
      +count(rownames(x[x[, "log2FoldChange"] > +SIG, ]) %in% goMap[, 1])
  }
  getDF2 <- function() {
    d[, "sums"] <<- rowSums(abs(d))
    d[, "group"] <<- factor(rownames(d), levels = rownames(d[order(d[, "sums"], decreasing = TRUE), ]))
    d[, "dntxt"] <<- -d[, "dn"]
    ## d[, "dnpc"] <- paste0(round(100 * (d[, "dntxt"] / length(sdeg$dn)), 1), "%")
    d[, "dnpc"] <<- paste0(round(100 * (d[, "dntxt"] / alldn), 1), "%")
    d[, "uptxt"] <<- +d[, "up"]
    ## d[, "uppc"] <- paste0(round(100 * (d[, "uptxt"] / length(sdeg$up)), 1), "%")
    d[, "uppc"] <<- paste0(round(100 * (d[, "uptxt"] / allup), 1), "%")
  }
  p <- function() {
    ggplot(d) +
      geom_col(aes(x = dn, y = group), fill = clrs[1]) +
      geom_col(aes(x = up, y = group), fill = clrs[2]) +
      geom_col(aes(x = dnannot, y = group), fill = clrs[4], width = 0.10, just = 0.5) +
      geom_col(aes(x = upannot, y = group), fill = clrs[4], width = 0.10, just = 0.5) +
      ## geom_text(aes(x = dn, y = group, label = dntxt), vjust = 0.5, hjust = -1) +
      ## geom_text(aes(x = up, y = group, label = uptxt), vjust = 0.5, hjust = +1) +
      geom_text(aes(x = dn, y = group, label = dntxt), vjust = 0.0, hjust = +1.1, size = 5) +
      geom_text(aes(x = up, y = group, label = uptxt), vjust = 0.0, hjust = -0.1, size = 5) +
      ## geom_text(aes(x = dn, y = group, label = dnannot), vjust = 0.5, hjust = +2.6, size = 4) +
      ## geom_text(aes(x = up, y = group, label = upannot), vjust = 0.5, hjust = -1.6, size = 4) +
      geom_text(aes(x = dn, y = group, label = dnpc), vjust = +1.5, hjust = +1.1, size = 5, color = "#444444") +
      geom_text(aes(x = up, y = group, label = uppc), vjust = +1.5, hjust = -0.1, size = 5, color = "#444444") +
      ggtitle(qq("@{t} SDEGs (total of @{alldn} down and @{allup} up, @{length(sdeg$names)} unique genes)")) +
      xlab("Counts and percent of all up/down-regulated SDEGs") +
      ylab("Group contrast (comparison vs. reference)") +
      xlim(-400, 400) +
      theme_std(plotDefs) +
      theme(
        legend.position = "NONE",
        text = element_text(size = 16),
        axis.text = element_text(size = 16, family = "Noto Sans Mono Condensed"),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 16, hjust = 1.0, vjust = 0.5))
  }

  alldn <- length(sdeg$dn)
  allup <- length(sdeg$up)
  clrs <- hcl.colors(4)
  n <- 1
  plotList <- list()
  for (t in toTitleCase(expFactors[["factorList"]][[expFactors[["factorList"]][["factors"]][1]]])) {
    d <- data.frame(dn = c(), up = c())
    lim <- 0
    for (i in lower1:upper1) {
      g1 <- expFactors[["factorList"]][["group"]][i]
      for (j in (lower2 + lim):upper2) {
        g2 <- expFactors[["factorList"]][["group"]][j]

        ## TODO: handle more than 2 factors at once
        for (f in 2:length(expFactors[["factorList"]][["factors"]])) {
          fvec <- expFactors[["factorList"]][[expFactors[["factorList"]][["factors"]][f]]]
          f1 <- fvec[which(paste0(fvec, t) %in% g1)]
          f2 <- fvec[which(paste0(fvec, t) %in% g2)]

          if(g1 %in% paste0(fvec, t) &&
             g2 %in% paste0(fvec, t)) {
            fmtstr <- "%16s"
            getDF1(sdegGroups[[i]][[j - 1 - lim]], f1, f2)
          }
        }
      }
      lim <- lim + 1
    }

    getDF2()
    plotList[[n]] <- p()
    n <- n + 1
  }

  d <- data.frame(dn = c(), up = c())
  for (k in 1:length(expFactors$factorList[[expFactors$factorList$factors[2]]])) {
    g1 <- paste0(expFactors$factorList[[expFactors$factorList$factors[2]]][k], toTitleCase(expFactors$factorList$time[2]))
    g2 <- paste0(expFactors$factorList[[expFactors$factorList$factors[2]]][k], toTitleCase(expFactors$factorList$time[1]))

    printf("Groups: %-20s %-20s\n", g1, g2)

    fmtstr <- "%20s"
    getDF1(sdegGroups[[i + 1]][[k]], f1 = g1, f2 = g2)
  }

  getDF2()
  plotList[[n]] <- (p() + ggtitle(qq("Late vs. Early SDEGs (total of @{alldn} down and @{allup} up, @{length(sdeg$names)} unique genes)")))

  return(plotList)
}

plotDiffExprAnnots <- function(sdegGroups, sdegAll) {
  for (t in expFactors$factorList$time) {
    d <- data.frame()
    t <- toTitleCase(t)
    lim <- 0
    for (i in lower1:upper1) {
      g1 <- expFactors[["factorList"]][["group"]][i]
      for (j in (lower2 + lim):upper2) {
        g2 <- expFactors[["factorList"]][["group"]][j]
        for (f in 2:length(expFactors$factorList$factors)) {
          fvec <- expFactors$factorList[[expFactors$factorList$factors[f]]]
          f1 <- fvec[which(paste0(fvec, t) %in% g1)]
          f2 <- fvec[which(paste0(fvec, t) %in% g2)]

          if(g1 %in% paste0(fvec, t) &&
             g2 %in% paste0(fvec, t)) {

            dn <- count(rownames(sdegGroups[[i]][[j - 1 - lim]][sdegGroups[[i]][[j - 1 - lim]][, "log2FoldChange" < -1], ]) %in% sdegAll)
            up <- count(rownames(sdegGroups[[i]][[j - 1 - lim]][sdegGroups[[i]][[j - 1 - lim]][, "log2FoldChange" > +1], ]) %in% sdegAll)

            fmtstr <- "%16s"
            d[paste0(sprintf(fmtstr, f1), sprintf(fmtstr, f2)), "dn"] <- dn
            d[paste0(sprintf(fmtstr, f1), sprintf(fmtstr, f2)), "up"] <- up
          }
        }
      }
      lim <- lim + 1
    }

    print(d)
  }
}

## plotTotalCounts(totalCountsAll, expFactors, plotDefs)
