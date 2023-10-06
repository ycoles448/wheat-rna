getArgs <- function() {
  return(commandArgs(trailingOnly = TRUE))
}

getConf <- function(file = "settings.toml") {
  return(read.config(file))
}
