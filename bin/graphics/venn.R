saveVennDiagrams <- function(prefix, plotList, legendList, expObjFactors, plotDefs) {
  d <- plotDefs[["defs"]]

  upList <- list()
  dnList <- list()
  for (i in 1:length(plotList)) {
    upList[[i]] <- plotList[[i]][[1]]
    dnList[[i]] <- plotList[[i]][[2]]
  }

  pVennUp <- arrangeGrob(
    do.call(arrangeGrob, c(upList, nrow = 2)),
    top = textGrob(
      paste(d["names.venn.sdeg"], "up-regulated", sep = ", "),
      gp = gpar(fontsize = d[["font.size.title"]])
    ),
    right = legendList[["up"]]
  )
  pVennDn <- arrangeGrob(
    do.call(arrangeGrob, c(dnList, nrow = 2)),
    top = textGrob(
      paste(d["names.venn.sdeg"], "down-regulated", sep = ", "),
      gp = gpar(fontsize = d[["font.size.title"]])
    ),
    right = legendList[["dn"]]
  )

  savePlotGrid(
    pVennUp,
    formatPath(prefix, "pVennUp.png", root = FALSE),
    plotDefs,
    width = dim(pVennUp[["grobs"]][[1]])[2] * 1500,
    height = dim(pVennUp[["grobs"]][[1]])[1] * 1500,
  )

  savePlotGrid(
    pVennDn,
    formatPath(prefix, "pVennDn.png", root = FALSE),
    plotDefs,
    width = dim(pVennDn[["grobs"]][[1]])[2] * 1500,
    height = dim(pVennDn[["grobs"]][[1]])[1] * 1500,
  )
}
