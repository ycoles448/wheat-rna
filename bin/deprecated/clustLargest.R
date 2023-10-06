# Only use Clust run with largest number of distinct modules
clust1 <- clustRuns[[which(counter == max(counter))]]
clust1[1:length(clust1) - 1, ] <- clust1[2:length(clust1), ]
clusters <- vector("list", ncol(clust1))
for (i in 1:ncol(clust1)) {
  colnames(clust1)[i] <- gsub("[.]genes[.]", "", gsub("C.*[.][.]", "", colnames(clust1)[i]))
  clusters[[i]] <- clust1[1:as.integer(colnames(clust1)[i]), i]
}

clust1list <- list()
clust1size <- c()
for (i in 1:ncol(clust1)) {
  clust1list[[i]] <- unlist(sapply(clust1[, i], function(x) {if (x != "") {x}}))
  names(clust1list[[i]]) <- NULL
  clust1size[i] <- length(clust1list[[i]])
}

clust1small <- clust1list[[which(clust1size == min(clust1size))]]
