formatPath <- function(..., root = TRUE) {
  # Argument checks
  if (!is.logical(root)) {
    stop("Argument root should be logical!")
  }

  argv <- c(as.list(environment()), list(...))
  argv <- unlist(argv)

  if (argv[1] == "TRUE" || argv[1] == "FALSE") {
    argv <- argv[-1]
  }

  ## if (OS == "Windows") {
  ## } else {

  # Check if path starts from the filesystem root
  if (root == TRUE) {
    if (substr(argv[[1]], 1, 1) == "/") {
      path <- character()
      for (i in argv) {
        path <- paste0(path, "/", i)
      }
      path <- substr(path, 2, nchar(path))
    } else {
      path <- getwd()
      for (i in argv) {
        path <- paste0(path, "/", i)
      }
    }
  } else {
    path <- character()
    for (i in argv) {
      path <- paste0(path, "/", i)
    }
    path <- substr(path, 2, nchar(path))
  }

  ## }

  return(path)
}

makePath <- function(path) {
  dir.create(path, recursive = TRUE)
}
