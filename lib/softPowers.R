exprMatr <- t(countsRL)

# Choose a set of soft-thresholding powers
powers <- c(seq(1, 10, 1))

# Call the network topology analysis function
sft <- pickSoftThreshold(exprMatr, powerVector = powers)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed rho-squared",
  type = "n",
  main = paste("Scale independence"),
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",

  ylab = "Mean Connectivity", type = "n",
  main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

dev.off()
