################################
################################
## Try permutation
################################
################################
require('Matrix')
library('foreach')
library('itertools')
library('doParallel')

nt <- detectCores()
cl <- makeCluster(nt)
registerDoParallel(cl)

results <- foreach(i = 1:nperm, .export = c('MD2S', 'ginv', 'irlba')) %dopar% {
  X1.perm <- sample.mat(sample.mat(X1)); y1.perm <- sample.mat(sample.mat(y1))

  b1.perm <- MD2S(X = X1.perm, y = y1.perm, dim = 5)
  out <- list("lz" = t(as.matrix(b1.perm$lz)),
              "lz.X" = t(as.matrix(b1.perm$lz.X)),
              "lz.y" = t(as.matrix(b1.perm$lz.y))
  )
  out
}
