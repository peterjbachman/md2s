# Purely internal functions, for use within md2s and md2sPermute

##standardize what's pased through, that's it.
## Give input
my.norm <- function(x) {
  #As vec
  x <- as.vector(x)
  #Norm each entry to mean
  x <- x - mean(x)
  #standardize by dividing above difference by the square root of the sum of entries squared
  (x / sum(x^2)^.5)
}

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

alpha.func <- function(x, z1, z2) {
  p1 <- exp(x) / (1 + exp(x))
  check.cor(p1 * z1 + (1 - p1) * z2)
}

alpha.func.X <- function(x, z1, z2) {
  p1 <- exp(x) / (1 + exp(x))
  check.cor.X(p1 * z1 + (1 - p1) * z2)
}

alpha.func.y <- function(x, z1, z2) {
  p1 <- exp(x) / (1 + exp(x))
  check.cor.y(p1 * z1 + (1 - p1) * z2)
}

alpha.func1 <- function(x) alpha.func(x, z.freq1, z.freq2)
alpha.func2 <- function(x) alpha.func(x, z.fit3, z.res3)
alpha.func2.X <- function(x) alpha.func.X(x, z.fit3, z.res3)
alpha.func2.y <- function(x) alpha.func.y(x, z.fit3, z.res3)

check.cor <- function(z.run) {
  z.run <- my.norm(z.run - mean(z.run))
  if (nrow(X) < ncol(X)) wX <- as.vector(XXprime %*% z.run) else wX <- as.vector(X1 %*% (t(X1) %*% z.run))
  if (nrow(y) < ncol(y)) wy <- as.vector(yyprime %*% z.run) else wy <- as.vector(y1 %*% (t(y1) %*% z.run))
  cov(wX, wy)
}

check.cor.X <- function(z.run) {
  z.run <- my.norm(z.run - mean(z.run))
  if (nrow(X) < ncol(X)) wX <- as.vector(XXprime %*% z.run) else wX <- as.vector(X %*% (t(X) %*% z.run))
  stats::var(wX)
}

check.cor.y <- function(z.run) {
  z.run <- my.norm(z.run - mean(z.run))
  if (nrow(y) < ncol(y)) wy <- as.vector(yyprime %*% z.run) else wy <- as.vector(y %*% (t(y) %*% z.run))
  stats::var(wy)
}

## Checks the correlation and outputs it?
## TODO Add inputs for x1 and y1 for my sanity
check.cor <- function(z.run) {
  wX <- as.vector(((X1) %*% t(X1)) %*% z.run)
  wy <- as.vector(((y1) %*% t(y1)) %*% z.run)
  stats::cov(wX, wy)
}

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
  # If there are fewer rows than columns
  if (nrow(z) <= ncol(z)) fits <- z %*% MASS::ginv(t(z) %*% z) %*% (t(z) %*% x)
  # If there are more columns than rows
  if (nrow(z) > ncol(z)) fits <- z %*% (MASS::ginv(t(z) %*% z) %*% t(z) %*% x)

  # Create residuals
  res <- x - fits

  # double-center the residuals (so the column, row, & grand means equal zero)
  res - make.int(res)
}

## Creates a sample matrix?
sample.mat <- function(X) {
  apply(X, 2, FUN = function(x) sample(x, length(x), replace = FALSE))
}

## No clue yet
tfidf <- function(mat) {
  tf <- mat

  ## show the number of columns that don't equal 0?
  id <- function(col) {
    sum(!col == 0)
  }

  idf <- log(nrow(mat) / apply(mat, 2, id))

  tfidf <- mat

  ## For every word
  for (word in names(idf)) {
    tfidf[, word] <- tf[, word] * idf[word]
  }
  return(tfidf)
}

## Normalizing function?
# `my.norm` creates z-scores for every value of input vector. - Patrick
#I'm not sure that it does. It's dividing by the square root of the sum of entries squared. It's just dividing by a common value that surpresses spread.
#I think that if this were creating a z-value, then we would have to multiply by root(n) and divide by some sd parameter.
my.norm <- function(x) {
  x <- as.vector(x)
  x <- x - mean(x)
  (x / sum(x^2)^.5)
}

##  impute missing values by iteratively double-centering the matrix. All na values are set to 0.
#  Note that the method choice for imputing values need not be of this method: Any method can be used.
dubcent.impute <- function(X) {
  X <- as.matrix(X)
  miss.mat <- is.na(X)
  X[is.na(X)] <- .5
  for (i in 1:10000) {
    Xlast <- X
    X <- X - make.int(X)
    X[miss.mat] <- 0
    if (sum(abs(Xlast - X)) == 0) break
  }
  X
}

## No clue yet
make.norm <- function(x) {
  edf <- sapply(x, FUN = function(z) sum(x <= z)) / (length(x) + 1)
  stats::qnorm(edf)
}
