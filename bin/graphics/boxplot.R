plotTotalCounts <- function(counts, expObjFactors, plotDefsObj) {
  d <- reshape2::melt(counts)
  defs <- plotDefsObj[["defs"]]
  colnames(d) <- "counts"

  d[, "group"] <- factor(
    expObjFactors[["factorList"]][["group"]],
    levels = expObjFactors[["factorList"]][["groupOrdered"]]
  )

  p <- ggplot(d, aes(x = group, y = counts, fill = group)) +
    geom_boxplot() +
    ggtitle(defs["names.barplot.totalCounts.title"]) +
    xlab(defs["names.xlab.cultivar"]) +
    ylab(defs["names.ylab.totalCounts"]) +
    theme_std(plotDefsObj) +
    guides_none()

  return(p)
}

## plotTotalCounts(totalCountsAll, expFactors, plotDefs)
