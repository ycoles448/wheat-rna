plotPairCountsScatter <- function(x, y, expObjFactors, plotDefsObj) {
  defs <- plotDefsObj[["defs"]]

  d <- data.frame(
    x = x,
    y = y,
    group = expObjFactors[["factorList"]][["group"]]
  )

  p <- ggplot(d, aes(x = x, y = y, color = group)) +
    geom_point() +
    stat_poly_line(se = TRUE, inherit.aes = FALSE, aes(x = x, y = y), color = "#00000040") +
    stat_poly_eq(vstep = 0.03, inherit.aes = FALSE, aes(x = x, y = y), color = "#00000080") +
    ggtitle(paste0(defs["names.scatter.pairCounts.title"], " (", expObjFactors[["factorList"]][["hkg"]], ")")) +
    xlab(paste0(defs["names.xlab.geneCounts"], ": ", expObjFactors[["factorList"]][["hkg"]])) +
    ylab(paste0(defs["names.ylab.totalCounts"])) +
    theme_std(plotDefsObj) +
    guides_none()

  return(p)
}

## plotTotalCounts(totalCountsAll, expFactors, plotDefs)
