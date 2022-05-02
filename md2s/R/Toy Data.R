# Toy data for 00_MD2S.R functions.

install.packages("matrixStats")
library(matrixStats)

# Input datasets
X <- matrix(rnorm(100*4, 30, 23), nrow = 100, ncol = 4)
y <- matrix(rnorm(100*8, 22, 5), nrow = 100, ncol = 8)

# create X1 and y1 in main file.

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
colSds(make.int(X))

#   Check row sd:
head(rowSds(X))
head(rowSds(make.int(X)))

#   Overall matrix sd:
sd(X)
sd(make.int(X))
