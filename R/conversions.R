
convertToColors <- function(mat, rng) {
  # Produce 'normalized' version of matrix, with values ranging from 0 to 1
  #rng <- range(mat, na.rm = TRUE)
  mat[mat<rng[1]] <- rng[1]
  mat[mat>rng[2]] <- rng[2]
  m <- (mat - rng[1])/diff(rng)
  m[is.nan(m)] <- NA
  # Convert to a matrix of sRGB color strings
  m2 <- m; class(m2) <- "character"
  #    m2[!is.na(m2)] <- rgb(colorRamp(rainbow(256))(m[!is.na(m)]), max = 255)
  m2[!is.na(m2)] <- rgb(colorRamp(c("white","midnightblue"))(m[!is.na(m)]), max = 255)
  m2[is.na(m2)] <- "transparent"
  return(m2)
}
