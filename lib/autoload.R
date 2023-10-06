#!/usr/bin/env Rscript

# Autoload.R, not intended to run directly, source from R scripts

# Author: Yutathkarn Coles
# Description: Script for automatically loading packages in R

# Does not take external configurations, requires a functional R installation
# path and appropriate compilers to install R packages

# Change to preferred repository URL
repo <- "https://cran.curtin.edu.au"

# Load dependencies from file
libs <- read.csv("rdeps.txt", sep = "\n", header = FALSE)

# Remove comments from library buffer
x <- c()
for (i in 1:length(libs[[1]])) {
  if (substr(libs[i, 1], 1, 1) == "#") {
    x <- append(x, i)
  }
}
libs <- libs[-x, 1] # Return vector from list
remove(x)

# Automatic package management
pkgLoad <- function(pkg) {
  ret <- tryCatch(
    {
      library(pkg, character.only = TRUE)
      return(0)
    },
    error = function(msg) {
      sprintf(as.character(msg))
      return(1)
    }
  )
  return(ret)
}

pkgInstall <- function(pkg, useBioconductor = TRUE) {
  ret <- tryCatch(
    {
      if (useBioconductor == TRUE) {
        invisible(BiocManager::install(pkgs = pkg, ask = FALSE))
      } else if (useBioconductor == FALSE) {
        install.packages(pkgs = pkg, repos = repo, quiet = TRUE)
      } else {
        sprintf("useBioconductor must be TRUE or FALSE\n")
        return(1)
      }
      return(0)
    },
    error = function(msg) {
      sprintf(as.character(msg))
      return(1)
    }
  )
  return(ret)
}

pkgSetupParallel <- function() {
  ret <- tryCatch(
    {
      pkgLoad("parallel")
      options(Ncpus = detectCores())
      return(0)
    },
    error = function(msg) {
      sprintf("Unable to load load parallel library\n")
      return(1)
    }
  )
  return(ret)
}

pkgSetupBiocManager <- function() {
  ret <- tryCatch(
    {
      pkgInstall("BiocManager", useBioconductor = FALSE)
      pkgLoad("BiocManager")
      return(0)
    },
    error = function(msg) {
      sprintf("Unable to install BiocManager\n")
      return(1)
    }
  )
  return(ret)
}

pkgSetupParallel()
pkgSetupBiocManager()

for (i in libs) {
  sprintf(paste0("Loading package: ", i, "\n"))
  if (pkgLoad(i) != 0) {
    sprintf(paste0("\nInstalling package: ", i, "\n"))
    if (pkgInstall(i) != 0) {
      sprintf(paste0("Could not install: ", i, "\n"))
    }
    pkgLoad(i)
  }
}
