getThreads <- function() {
  tryCatch(
    {
      threads <- conf[["general"]][["threads"]]
      if (threads == 0) {
        threads <- detectCores()
      }
    },
    error = function() {
      message("Could not get threads from config")
      threads <- detectCores()
    }
  )

  return(threads)
}

threads <- getThreads()
registerDoParallel(cores = threads)
register(MulticoreParam(workers = threads))

if ("WGCNA" %in% installed.packages()) enableWGCNAThreads(nThreads = threads)
