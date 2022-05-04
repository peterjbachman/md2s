# Toy data for 00_MD2S.R functions.

library(matrixStats)

# Input datasets
# create X1 and y1 in main file.
X <- matrix(rnorm(100*4, 30, 23), nrow = 100, ncol = 4)
y <- matrix(rnorm(100*8, 22, 5), nrow = 100, ncol = 8)


# WHAT IS `make.int` FUNCTION DOING???
#   SOLUTION: make.int() standardizes the columns & rows of the matrix:
#   "we preprocess the matrices by double-centering them, so that the row-mean,
#   column-mean, and grand mean is zero" (pg. 216 of paper)

#   Check column means:
colMeans(X)
colMeans(make.int(X))
#   Check row means:
head(rowMeans(X))
head(rowMeans(make.int(X)))
#   Overall matrix mean:
mean(X)
mean(make.int(X))

#   Check column sd:
colSds(X)
colSds(make.int(X)) # Now all columns have same std.
#   Check row sd:
head(rowSds(X))
head(rowSds(make.int(X))) # Now all rows have same std.
#   Overall matrix sd:
sd(X)
sd(make.int(X))

# Now check if column/row/grand means are zero:
colMeans(X1) # Basically zero.
head(rowMeans(X1)) # basically zero.
mean(X1) # Basically zero.
colMeans(y1) # Basically zero.
head(rowMeans(y1)) # basically zero.
mean(y1) # Basically zero.



# WHAT IS `my.norm` FUNCTION DOING?
#   ANSWER:  `my.norm` creates z-scores for every value of input.
test1 <- t(b0$z%*%X1)
my.norm <- function(x) {
  x <- as.vector(x)
  x <- x - mean(x)
  (x / sum(x^2)^.5)
}
# Test with
