#' Permutation Test
#'
#' Permutation Test
#'
#'
#' @param kX.num K(2) as defined in the paper.
#' @param n Number of rows in the matrix
#' @param ky K(1) in the paper
#' @param nsims The number of simulations to be done. Default is set to 200.
#' @param nperm The number of permutations to be done. Default is set to 200.
#' @param nboot The number of bootstraps to be performed. Default is set to 200.
#' @param d number of dimensions
#' @param fpath Designated file path where we will store the results
#'
#' @return An object containing:
#' \item{kX.num}{Number of columns in the matrix, K(2) in the paper}
#' \item{n}{Number of rows in the matrix}
#' \item{ky}{Number of columns of matrix, K(1) in the paper}
#' \item{nsims}{Number of simulations done}
#' \item{nperm}{Number of permutations done}
#' \item{nboot}{Number of bootstraps done}
#' \item{fpath}{Designated file path for the result}
#'
#' @author Cecilia Y. Sui and Evan E. Jo: \email{c.sui@@wustl.edu} \email{ejo@wustl.edu}
#' @seealso [md2s::md2s()]
#' @rdname md2sPermute
#' @include md2sPermute.R
#' @aliases md2spermute
#' @importFrom foreach %dopar%
#'
#' @examples
#' \dontrun{
#' md2sPermute(kX.num = 100, n = 50, ky = 40, nsims = 200, nperm = 200, nboot = 200)
#' }
#'
#'
#'
#'
#'
#' @export
md2sPermute <- function(kX.num = 100, n = 50, ky = 40, nsims, nperm, nboot, fpath = ".", sim = TRUE) {
  ## Generate Data
  results.all <- NULL
  for (j in 1:nsims) {
    for (kX in c(kX.num)) {

      if (sim) {
      ## -------------------------
      ## Create z.s y z2.s
      ## -------------------------
        z1 <- stats::rnorm(n)
        z2 <- stats::rnorm(n)
        zs.true <- cbind(z1, z2)
        zX.true <- cbind(stats::rnorm(n), stats::rnorm(n), stats::rnorm(n))
        zy.true <- cbind(stats::rnorm(n), stats::rnorm(n))

        bsX.true <- cbind(stats::rnorm(kX), stats::rnorm(kX))
        bX.true <- cbind(stats::rnorm(kX), stats::rnorm(kX), stats::rnorm(kX))
        bsy.true <- cbind(stats::rnorm(ky), stats::rnorm(ky))
        by.true <- cbind(stats::rnorm(ky), stats::rnorm(ky))

        dw.z <- diag(c(2, 1))
        dw.y <- diag(c(4, 2))
        dw.X <- diag(c(4, 2, 2))

        X1 <- zs.true %*% dw.z %*% t(bsX.true) + zX.true %*% dw.X %*% t(bX.true) + 2 * matrix(stats::rnorm(n * kX), nrow = n)
        y1 <- zs.true %*% dw.z %*% t(bsy.true) + zy.true %*% dw.y %*% t(by.true) + 2 * matrix(stats::rnorm(n * ky), nrow = n)
      }


      # Notice here b1 is a MD2S obejct #
      b1 <- md2s(X = (X1), y = (y1), dim = 5)

      ################################
      ## Try permutation
      ################################

      nt <- parallel::detectCores()
      cl <- parallel::makeCluster(nt)

      # delete ginv if getting errors on this
      parallel::clusterExport(cl, c("md2s", "irlba", "sample.mat",
                                    "make.int", "fastres", "md2sInner",
                                    "cleanup", "my.norm", "check.cor", "alpha.func"))
      doParallel::registerDoParallel(cl)
      out <- Matrix::Matrix()

      results <- foreach::foreach(i = 1:nperm, .packages = "MASS", "parallel") %dopar% {
        X1.perm <- sample.mat(sample.mat(X1))
        y1.perm <- sample.mat(sample.mat(y1))
        b1.perm <- md2s(X = X1.perm, y = y1.perm, dim = 5)
        out <- list(
          "lz" = t(as.matrix(b1.perm$lz)),
          "lz.X" = t(as.matrix(b1.perm$lz.X)),
          "lz.y" = t(as.matrix(b1.perm$lz.y))
        )
        out
      }

      foreach::registerDoSEQ()
      parallel::stopCluster(cl)

      lz.run <- lz.X.run <- lz.y.run <- list()
      for (i in 1:length(results)) {
        lz.run[[i]] <- results[[i]]$lz
        lz.X.run[[i]] <- results[[i]]$lz.X
        lz.y.run[[i]] <- results[[i]]$lz.y
      }

      lz.run2 <- do.call("rbind", lz.run)
      lz.X.run2 <- do.call("rbind", lz.X.run)
      lz.y.run2 <- do.call("rbind", lz.y.run)

      pz1 <- pz2 <- pz2.y <- pz2.X <- NULL
      dims <- 5

      for (i in 1:dims) {
        pz2[i] <- mean(b1$lz[i] < lz.run2[, i])
        pz2.y[i] <- mean(b1$lz.y[i] < lz.y.run2[, i])
        pz2.X[i] <- mean(b1$lz.X[i] < lz.X.run2[, i])
      }

      pvals.s <- pz2
      pvals.X <- pz2.X
      pvals.y <- pz2.y

      scenario2 <- 1
      results.curr <- c(n, kX, scenario2, pvals.s, pvals.X, pvals.y)

      names(results.curr) <- c("N", "kX", "scen2", paste("ps_", 1:5, sep = ""), paste("pX_", 1:5, sep = ""), paste("py_", 1:5, sep = ""))

      results.all <- results.curr

      # return(results.all)
      ran.name <- paste0(fpath, "output_", round(stats::runif(1), 10) * 1e10, "_permtest.RData")
      save(results.all, file = ran.name)
    }
  }


  ## Aggregate results
  files.all <- list.files(fpath, pattern = "output_")
  if(!exists("out.all")){
    out.all<-NULL
    for(i in files.all) {
      load(paste(fpath, i, sep = ""))
      out.all<-rbind(out.all,results.all)
    }
  }
  save(out.all, file = paste0(fpath, "outall.RData"))

  for(i in files.all) {
    file.remove(paste(fpath, i, sep = ""))
  }

}

