expObj <- setRefClass(
  "expObj",
  fields = list(
    mat = "data.frame",
    meta = "data.frame",
    contrast = "character",
    factors = "vector",
    nfactors = "integer",
    nreps = "integer",
    species = "character"
  )
)

expObjFactors <- setRefClass(
  "expObjFactors",
  fields = list(
    factorList = "list"
  )
)

plotDefsObj <- setRefClass(
  "plotDefsObj",
  fields = list(
    defs = "vector"
  )
)

geneListObj <- setRefClass(
  "geneListObj",
  fields = list(
    genes = "vector"
  )
)
