getContrastsFungi <- function(dds, expObj = expObj) {
  pairs <- list()

  for (i in 1:length(cultivars)) {
    g1 <- paste0(toTitleCase(cultivars[i]), toTitleCase(times[2]))
    g2 <- paste0(toTitleCase(cultivars[i]), toTitleCase(times[1]))

    res.dds <- results(
      dds,
      alpha = 0.05,
      contrast = c(CONTRAST, g1, g2)
    )

    pairs[[i]] <- res.dds
  }

  return(pairs)
}
