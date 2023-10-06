# Defining standards for plots
getPlotDefs <- function(defs) {
  return(plotDefsObj(defs = unlist(defs)))
}

theme_std <- function(plotDefsObj) {
  d <- plotDefsObj[["defs"]]
  t <- theme_classic() +
    theme(
      axis.text.x = element_text(
        size = as.numeric(d["font.size.axis"]),
        angle = as.numeric(d["font.angle.x"]),
        hjust = 1,
      ),
      axis.text.y = element_text(
        size = as.numeric(d["font.size.axis"]),
        angle = as.numeric(d["font.angle.y"]),
      )
    )

  return(t)
}

theme_empty <- function(plotDefsObj) {
  d <- plotDefsObj[["defs"]]
  t <- theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
  )
  return(t)
}

theme_venn <- function(plotDefsObj) {
  d <- plotDefsObj[["defs"]]
  m <- d["venn.margin"]

  t <- theme(
    plot.margin = margin(m, m, m, m, unit = "cm")
  )
  return(t)
}

guides_none <- function() {
  t <- theme(
    legend.position = "none"
  )

  return(t)
}

savePlotStd <- function(plot, filename, plotDefsObj, outDir = conf[["data"]][["plots"]], width = NA, height = NA) {
  file <- formatPath(outDir, filename)
  makePath(dirname(file))

  d <- plotDefsObj[["defs"]]

  if (is.na(width)) {
    x <- as.numeric(d["size.std.x"]) * cols
  } else {
    x <- as.numeric(width)
  }

  if (is.na(height)) {
    y <- as.numeric(d["size.std.y"]) * rows
  } else {
    y <- as.numeric(height)
  }
  u <- "px"

  print(paste0("> Saving plot to ", file, " with dimensions ", x, ", ", y))

  ggsave(
    filename = file,
    plot = plot,
    units = u,
    width = x,
    height = y,
  )
}

savePlotGrid <- function(plot, filename, plotDefsObj, outDir = conf[["data"]][["plots"]], width = NA, height = NA) {
  file <- formatPath(outDir, filename)
  makePath(dirname(file))

  d <- plotDefsObj[["defs"]]

  rows <- dim(plot)[1]
  cols <- dim(plot)[2]

  if (is.na(width)) {
    x <- as.numeric(d["size.std.x"]) * cols
  } else {
    x <- as.numeric(width)
  }

  if (is.na(height)) {
    y <- as.numeric(d["size.std.y"]) * rows
  } else {
    y <- as.numeric(height)
  }
  u <- "px"

  print(paste0("> Saving plot to ", file, " with dimensions ", x, ", ", y))

  ggsave(
    filename = file,
    plot = plot,
    units = u,
    width = x,
    height = y,
    limitsize = FALSE,
  )
}

saveHeatmap <- function(plot, filename, plotDefsObj, outDir = conf[["data"]][["plots"]], width = NA, height = NA) {
  file <- formatPath(outDir, filename)
  makePath(dirname(file))

  d <- plotDefsObj[["defs"]]

  if (is.na(width)) {
    x <- as.numeric(dim(attr(p3, "matrix"))[1]) * 10 + 200
  } else {
    x <- width
  }

  if (is.na(height)) {
    y <- as.numeric(dim(attr(p3, "matrix"))[1]) * 10 + 200
  } else {
    y <- height
  }
  u <- "px"

  printf("x: %s\ty: %s\n", x, y)

  png(
    filename = file,
    units = u,
    width = x,
    height = y,
  )

  draw(plot)
  dev.off()
}
