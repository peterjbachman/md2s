# Purely internal functions, for use within md2s and md2sPermute

# Used in md2sInner, md2s,
##standardize what's passed through, that's it.
## Give input
my.norm <- function(x) {
  #As vec
  x <- as.vector(x)
  #Norm each entry to mean
  x <- x - mean(x)
  #standardize by dividing above difference by the square root of the sum of entries squared
  (x / sum(x^2)^.5)
}

# Only used in md2sInner
##Input optional covariates for explaining scaled locations in the shared subspace
cleanup <- function(X.c) {
  #If dim of X.c is greater than 0, then...
  if (length(X.c) > 0) {
    #apply my.norm then,
    X.c <- apply(X.c, 2, FUN = function(x) x - mean(x))
    #get square of col means greater than below number. Then, store as X.c
    X.c <- X.c[, colMeans(X.c^2) > 1e-10]
  }
  X.c
}

alpha.func <- function(x, z1, z2, option, XXprime, X1, yyprime, y1) {
  p1 <- exp(x) / (1 + exp(x))
  check.cor(p1 * z1 + (1 - p1) * z2, option = option, XXprime = XXprime, X1 = X1, yyprime = yyprime, y1 = y1)
}

check.cor <- function(z.run, option, XXprime, X1, yyprime, y1) {
  z.run <- my.norm(z.run - mean(z.run))
  if (nrow(X1) < ncol(X1)) wX <- as.vector(XXprime %*% z.run) else wX <- as.vector(X1 %*% (t(X1) %*% z.run))
  if (nrow(y1) < ncol(y1)) wy <- as.vector(yyprime %*% z.run) else wy <- as.vector(y1 %*% (t(y1) %*% z.run))

  if (option == "X") {
    stats::var(wX)
  } else if (option == "y") {
    stats::var(wy)
  } else {
    stats::cov(wX, wy)
  }
}

# Used in md2s
# `make.int` standardizes the rows (columns) of a matrix such that each
#   row (column) of the matrix has the same standard deviation & mean. I think
#   it also does something with the grand mean.
make.int <- function(X) {
  int1 <- rep(1, nrow(X)) %*% t(colMeans(X))
  int2 <- rowMeans(X) %*% t(rep(1, ncol(X)))
  int1 + int2 - mean(int1 + int2) + mean(X)
}

## Create residuals faster?
fastres <- function(x, z) {
  z <- cbind(1, z)
  #   NOTE: these following steps create fitted values.
  #   If there are fewer rows than columns
  if (nrow(z) <= ncol(z)) fits <- z %*% MASS::ginv(t(z) %*% z) %*% (t(z) %*% x)
  # If there are more columns than rows
  if (nrow(z) > ncol(z)) fits <- z %*% (MASS::ginv(t(z) %*% z) %*% t(z) %*% x)

  # Create residuals
  res <- x - fits

  # double-center the residuals (so the column, row, & grand means equal zero)
  res - make.int(res)
}

#
## Creates a sample matrix?
sample.mat <- function(X) {
  apply(X, 2, FUN = function(x) sample(x, length(x), replace = FALSE))
}
