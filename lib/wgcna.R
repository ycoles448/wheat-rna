## # WGCNA
## ## exportTraitsFungi(formatPath(conf[["data"]][["data"]], paste0(expData[["species"]], "Traits.tsv"), root = FALSE), meta[meta[, "control"] == "disease", ], c(factors, extraFactors))
## exportTraitsFungi(formatPath(conf[["data"]][["data"]], paste0(expData[["species"]], "Traits.tsv"), root = FALSE), meta[meta[, "control"] == "disease" & meta[, "time"] == "late", ], c(factors, extraFactors))

## ## exprData <- t(countsRL)
## exprData <- t(countsRL[, rownames(meta[meta[, "control"] == "disease" & meta[, "time"] == "late", ])])
## ## exprData <- fixDataStructure(exprData)
## traits <- read.csv(fileTraits, header = TRUE, sep = "\t")
## if (colnames(traits)[length(traits)] == "X") {
##   traits <- traits[, -c(length(traits))]
## }
## traits <- traits[, -1]
## rownames(traits) <- rownames(exprData)

## t <- traits[, c("rating.SVS", "rating.MS", "rating.MR", "rating.R")]

## # TODO: Implement proper checkpointing
## ## x <- list.files(formatPath(conf[["data"]][["data"]], conf[["data"]][["tom"]]))
## ## if (length(x) > 0) {
## ##   for (i in x) load(formatPath(conf[["data"]][["data"]], conf[["data"]][["tom"]], i))
## ## } else {
## cor <- WGCNA::cor

## # See below for manual heirachial clustering, blockwiseModules (and similar) recommended for large datasets
## exprData <- fixDataStructure(exprData)
## net <- blockwiseConsensusModules(
##   exprData,
##   maxBlockSize = 20000,
##   power = 6,
##   minModulesSize = 30,
##   deepSplit = 2,
##   pamRespectsDendro = TRUE,
##   mergeCutHeight = 0.25,
##   numericLabels = TRUE,
##   minKMEtoStay = 0,
##   saveTOMs = FALSE,
##   verbose = 5,
## )

## ## col = labels2colors(net[["colors"]])
## ## plotDendroAndColors(
## ##   net[["dendrograms"]][[1]],
## ##   col[net[["blockGenes"]][[1]]],
## ##   "Module colours",
## ##   groupLabels = "Gene networks",
## ##   dendroLabels = FALSE,
## ##   hang = 0.03,
## ##   addGuide = TRUE,
## ##   guideHang = 0.05,
## ## )

## ## sizeGrWindow(9, 5)

## ## TODO: Take unaligned reads from wheat and ptr, for Leon's metagenomics project

## ## Make plots of reads that aligned and did not align
