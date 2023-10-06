clr2rgb <- function(col) {
  x <- col
  r <- sqrt(strtoi(paste0("0x", substr(x, 2, 3))) / 255)
  g <- sqrt(strtoi(paste0("0x", substr(x, 4, 5))) / 255)
  b <- sqrt(strtoi(paste0("0x", substr(x, 6, 7))) / 255)

  return(c(r, g, b))
}
