# Configs
argv <- getArgs()
conf <- getConf()

# Fungi analysis
SIG <- 0.01
L2FC <- 2
countThreshold <- 50

species <- "ptr"
confPlots <- getConf("lib/graphics.toml")

fileMat <- conf[["fungi"]][["counts"]]
fileMeta <- conf[["fungi"]][["meta"]]
fileTraits <- conf[["fungi"]][["traits"]]
fileGoMap <- formatPath(conf[["data"]][["data"]], conf[["data"]][["go"]], paste(species, "map-go.txt", sep = "-"))
# fileGoNames is dynamically loaded

dirClust <- formatPath(conf[["data"]][["data"]], conf[["data"]][["clust"]], species)

mat <- read.csv(fileMat, sep = "\t", header = TRUE)
meta <- read.csv(fileMeta, sep = "\t", header = TRUE)

genes <- conf[["data"]][["genesFungi"]]
genesSelect <- unlist(getConf("genes.toml")[["fungi"]])
hkg <- conf[["fungi"]][["hkg"]]

contrast <- conf[["analyses"]][["contrast"]]
factors <- conf[["fungi"]][["factors"]]
extraFactors <- conf[["fungi"]][["extraFactors"]]
factorList <- conf[["factors"]]

nreps <- conf[["fungi"]][["nreps"]]
samples <- which(meta[["control"]] != "control")

clustRuns <- vector("list", length(list.files(dirClust)))
names(clustRuns) <- list.files(dirClust)
ont <- c("BP", "MF", "CC")
