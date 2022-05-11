#' md2s
#'
#' Multi-Dataset Multidimensional Scaling
#'
#'
#' @param X First matrix of realized outcomes \eqn{m = 1} for \eqn{N} observations and \eqn{K_1} covariates; dimensions \eqn{N \times K_1}.
#' @param y Second matrix of realized outcomes \eqn{m = 2} for \eqn{N} observations and \eqn{K_2} covariates; dimensions \eqn{N \times K_1}
#' @param X.s Optional matrix of covariates for estimating scaled locations in the shared subspace.
#' @param X.X Optional matrix of covariates for estimating scaled locations in the idiosyncratic subspace of `X`.
#' @param X.y Optional matrix of covariates for estimating scaled locations in the idiosyncratic subspace of `y`.
#' @param init INOPERABLE. Matrix initialization method.
#' @param dim Number of dimensions to be estimated by the model (from most- to least-explanatory).
#' @param tol Iteration tolerance threshold. Criteria for termination of iteration.
#' @param BIC INOPERABLE. Bayesian Information Criterion.
#'
#' @return An object containing:
#' \item{z}{TODO}
#' \item{z.X}{TODO}
#' \item{z.y}{TODO}
#' \item{w.Xs}{TODO}
#' \item{w.ys}{TODO}
#' \item{beta.z}{TODO}
#' \item{proportionX}{TODO}
#'
#' @author Peter Bachman <bachman.p@wustl.edu>, Patrick Edwards <edwards.p@wustl.edu>, and Zion Little <l.zion@wustl.edu>
#' @seealso [md2s::md2sPermute()]
#' @include md2s.R
#'
#' @examples
#' \dontrun{
#' md2s()
#' }
#'
#'
#'
#'
#
#' @export

md2s <- function( # List of arguments and descriptions below:
                 X, # $N \mult K_1$ dataset $m = 1$ of realized outcomes.
                 #  DESCRIPTION: first dataset of realized outcomes for `N` observations and `K_1` covariates.
                 #    Note: $K_1$ & $K_2$ can differ.
                 #  EXAMPLE: roll-call votes for `N` legislators with `K_1` covariates.

                 y, # $N \mult K_2$ dataset $m = 2$ of realized outcomes.
                 #  DESCRIPTION: second dataset of realized outcomes for `N` observations and `K_2` covariates.
                 #    Note: `K_1` & `K_2` can differ.
                 #  EXAMPLE: floor speech text for `N` legislators with `K_2` covariates.

                 X.s = NULL, # Optional covariates for estimating scaled locations in shared subspace.
                 #  DESCRIPTION: dataset of covariates that explain observations' scaled locations in shared subspace.
                 #  EXAMPLE: covariates that explain legislators' general ideological placements.

                 X.X = NULL, # Optional covariates for estimating scaled locations in idiosyncratic subspace of X.
                 #  DESCRIPTION: dataset of covariates that explain observations' scaled locations in idiosyncratic subspace of X.
                 #  EXAMPLE: covariates that explain legislators' ideological placements as determined by roll call votes.

                 X.y = NULL, # Optional covariates for estimating scaled locations in idiosyncratic subspace of y.
                 #  DESCRIPTION: dataset of covariates that explain observations' scaled locations in idiosyncratic subspace of y.
                 #  EXAMPLE: covariates that explain legislators' ideological placements as determined by floor speech text.

                 init = "svd", # INOPERABLE. Matrix initialization method.
                 #  NOTE 1: Default is Singular Value Decomposition (SVD) of a matrix.
                 #  NOTE 2: Function permits no other argument besides default "svd".

                 dim = 1, # Number of dimensions to be fitted by the model (from most- to least-explanatory).
                 #  [BETTER DESCRIPTION NEEDED]

                 tol = 1e-6, # TOL = convergence tolerance parameter.
                 #  DESCRIPTION: criteria for iteration termination.
                 #    Iterations terminate if difference between current & prior iteration are less than `tol` value.

                 BIC = FALSE # INOPERABLE. Bayesian Information Criterion (?)
) {

  # Sets `bic.sort` to `NULL` if the `BIC` argument equals `FALSE`, which sets the "bic" element of the `output` list to `NULL`
  if (BIC == FALSE) bic.sort <- NULL

  # Renames optional covariates for explaining scaled locations in the shared subspace.
  X.c <- X.s

  # Sets `n` as the number of observations in X, the first dataset (m = 1) of realized outcomes.
  n <- nrow(X)

  # Creates three `Z` column vectors of 1s with length = number of rows in X.
  #   Correspond to the shared (Z_S) & idiosyncratic (Z_{(M)}) subspaces.
  #     (from paper): Z_S contains latent locations in the shared subspace in columns for each `dim` dimensions (in this case, 1).
  #     (from paper): Z_{(M)} contains latent locations in the idiosyncratic subspace for `dim` latent dimensions (in this case, 1).
  ZX.mat <- Zy.mat <- Z.mat <- as.matrix(rep(1, n))


  # Creates one vector of 1s with length = $K_1$, the number of covariates in first matrix of realized outcomes $X$ (i.e., `ncol(X)`).
  #   (from paper): W_{(M)} is matrix of shared factors for the shares subspace for the dataset Y_{(M)}.
  wXs.mat <- wX.mat <- rep(1, ncol(X))

  # Creates one vector of 1s with length = $K_2$, the number of covariates in second matrix of realized outcomes $y$ (i.e., `ncol(y)`).
  #   (from paper): W_{(M)} is matrix of shared factors for the shared subspace for the dataset Y_{(M)}.
  wys.mat <- wy.mat <- rep(1, ncol(y))

  # REPEAT of third to last line of code.
  #   Delete unnecessary repeat?
  ZX.mat <- Zy.mat <- Z.mat <- as.matrix(rep(1, n))

  # Creates `NULL` objects (??for later use??).
  #   Potential purpose: proportionX may be intended to measure how much larger the `X` dataset is than the `y` dataset. It does nothing here, though.
  proportionX <- lz <- lz.X <- lz.y <- NULL

  # Double-center `X1` and `y1`:
  #   "we preprocess the matrices by double-centering them, so that the row-mean, column-mean, and grand mean is zero" (pg. 216 of paper).
  #   BACKGROUND: the `make.int` functio
  # `make.int` standardizes the rows (columns) of a matrix such that each row (column) of the matrix has the same standard deviation & mean.
  X1 <- X - make.int(X)
  y1 <- y - make.int(y)

  #
  for (i in 1:dim) {

    # Concatenates and prints "Fitting dimension `i`" for each `i`.
    cat("----  Fitting Dimension ", i, "----  \n")

    # FUNCTION: `fastres` function defined below.
    #   BACKGROUND: `fastres` creates double-centered residuals.
    X1 <- fastres(
      x = X1, # double-centered matrix for first input dataset (m = 1)
      z = cbind( # cbind for two same-length column vectors filled with 1s.
        Z.mat, # column vector (see above); latent factors of shared subspace.
        ZX.mat # column vector (see above); latent factors for subspace of 1st input dataset (m = 1)
      ) # Here, `fastres` is creating fitted values by multiplying \beta = 3 x N matrix with X1.
    ) # output is similar to original X1 but now the row/column/grand means are EVEN CLOSER to zero.

    # Same as above for y1:
    y1 <- fastres(y1, cbind(Z.mat, Zy.mat))

    # Same as above except performed on t(X1) and with column vectors of 1s with length = # of X covariates.
    X1 <- t(fastres(
      t(X1),
      cbind(
        wX.mat, # matrix of 1s of length ncol(X) (i.e., # of covariates in 1st dataset m = 1)
        wXs.mat # Same as wX.mat
      )
    )) # Again, the means are even closer to zero generally.

    # Same as directly above for y1:
    y1 <- t(fastres(t(y1), cbind(wy.mat, wys.mat)))

    # double-centers new X1:
    #   Again, means of new result is basically closer to zero.
    X1 <- X1 - make.int(X1)

    # Same as directly above but for new y1:
    y1 <- y1 - make.int(y1)

    # (???) WHAT DOES `MD2S_inner` FUNCTION DO (???)
    b0 <- md2sInner(
      X0 = X1, # Double-centered/scaled, derived from X.
      y0 = y1, # Double-centered/scaled, derived from y1.
      X.c.0 = X.c, # covariates associated with shared subspace.
      X.X.0 = X.X, # Covariates associated with X subspace.
      X.y.0 = X.y, # covariates associated with y subspace.
      tol0 = tol # tolerance parameter.
    )
    # ,burnin=0,gibbs=100,thin=1) - ORIGINAL COMMENT! Potentially additional arguments they planned to provide??

    Z.mat <- cbind(
      Z.mat, # (100 x 1) col. vector of 1s; (LIKELY) intercept term.
      b0$z # (100 x 1) col. vector; (LIKELY) latent locations on shared subspace.
      # DESCRIPTION: (from paper) "Z_S contains latent locations on the
      #   shared subspace in columns for each of the Q_S dimensions" (pg. 216).
      #   Likely only 1 column because `dim` = 1. If `dim` > 1 then `ncol(b0$z)` > 1.
      # EXAMPLE: vector of scaled ideological locations of legislators.
    ) # OUTPUT: (N x 2) matrix containing column of 1s & column of z's.

    Zy.mat <- cbind(
      Zy.mat, # (N x 1) col. vector of 1s; (LIKELY) intercept term.
      b0$z.y # (N x 1) col. vector; (LIKELY) latent locations on y's subspace.
      # DESCRIPTION: (from paper) "each column of Z_{(m)} contains the latent
      #   locations in the idiosyncratic subspace for Q_{(m)} latent dimensions" (pg. 216).
      #   Again, likely only 1 column because `dim` = 1.
      # EXAMPLE: vector of scaled ideological locations of legislators based
      #   on roll call votes.
    ) # OUTPUT: (N x 2) matrix containing column of 1s & column of z.y's.

    ZX.mat <- cbind(
      ZX.mat, # (100 x 1) col. vector of 1s; (LIKELY) intercept term.
      b0$z.X # (100 x 1) col. vector; (LIKELY) latent locations on X's subspace.
      # DESCRIPTION: (from paper) "each column of Z_{(m)} contains the latent
      #   locations in the idiosyncratic subspace for Q_{(m)} latent dimensions" (pg. 216).
      #   Again, likely only 1 column because `dim` = 1.
      # EXAMPLE: vector of scaled ideological locations of legislators based
      #   on floor speech text data.
    ) # OUTPUT: (N x 2) matrix containing column of 1s & column of z.X's


    # OVERVIEW: (LIKELY) the next four 'commands' of code normalize (or standardize
    #   to z-scores that follow the normal dist.) the latent locations on each relevant subspace.

    wXs.mat <- cbind(
      wXs.mat, # A (row?) vector of 1s with length = K_1 (# of covariates in X, or `ncol(X)`).
      my.norm( # `my.norm` creates z-scores for every value of input.
        t(b0$z %*% X1)
      ) # OUTPUT: (K_1 x 1) column vector.
    ) # OUTPUT: (K_1 x 2) matrix with a column of 1s and column of z-scores.
    # (from paper): "W_{(M)} is matrix of shared factors for the shared
    #   subspace for the dataset Y_{(M)}" (pg. 216).

    # Same as above except using `y1`.
    wys.mat <- cbind(wys.mat, my.norm(t(b0$z %*% y1)))

    # Same as above except using `z.X` and `X1`.
    wX.mat <- cbind(wX.mat, my.norm(t(b0$z.X %*% X1)))

    # Same as above except using `z.y` and `y1`.
    wy.mat <- cbind(wy.mat, my.norm(t(b0$z.y %*% y1)))


    # (LIKELY) this is to assess if two L_{(m)} are proportional:
    #   (from paper) "We assume that the two matrices L_{(1)} and L_{(2)} are
    #     proportional, so any difference between them is attributable to the
    #     relative scales across data sources Y_{(m)}
    proportionX[i] <- b0$pr # Recall `i in dim`, so `i` is a dimension.


    # (LIKELY) creates the square L_S matrix with `ncol` = `nrow` = `dim`,
    #   or the number of latent dimensions of the shared subspace Q_S.
    # (From paper) "L_S is Q_S ? Q_S non-negative, diagonal matrix of loadings
    #   for the shared subspace" (pg. 216).
    if (nrow(X1) < max(
      ncol(X1), # Equals K_1 (# of covariates in X)
      ncol(y1) # Equals K_2 (# of covariates in y)
    )
    ) {
      lz[i] <- t(b0$z) %*% (X1 %*% t(X1)) %*% (y1 %*% t(y1)) %*% b0$z
    } else {
      lz[i] <- ((t(b0$z) %*% X1) %*% t(X1)) %*% (y1 %*% (t(y1) %*% b0$z))
    } # OUTPUT: (1 x 1) matrix (for single dimension, `dim` = 1).

    # (LIKELY) following two lines create the square L_{(m)} matrix with `ncol` = `nrow` = `dim`,
    #   or the number of latent dimensions of the idiosyncratic subspace Q_{(m)}.
    # (From paper) "L_S is Q_{(m)} ? Q_{(m)} non-negative, diagonal matrix of
    #   loadings for the idiosyncratic subspaces" (pg. 216).
    lz.X[i] <- sum((t(X1) %*% b0$z.X)^2)
    lz.y[i] <- sum((t(y1) %*% b0$z.y)^2)
  } # For-loop repeats for every dimension estimated (i.e., when `dim` > 1)

  # What the heck? This is literally just `Z.mat[,2]`
  #   Check for yourself by running:
  #   `as.matrix(Z.mat[,2]) == as.matrix(Z.mat[,apply(Z.mat,2,sd) > 0])`
  #   Everything is `TRUE`. Why do it like this?
  Z.mat <- as.matrix(Z.mat[, apply(Z.mat, 2, stats::sd) > 0])

  # Same as above but for Zy.mat:
  Zy.mat <- as.matrix(Zy.mat[, apply(Zy.mat, 2, stats::sd) > 0])

  # Same as above but for ZX.mat:
  ZX.mat <- as.matrix(ZX.mat[, apply(ZX.mat, 2, stats::sd) > 0])

  # Same as above but for wX.mat:
  wX.mat <- as.matrix(wX.mat[, apply(wX.mat, 2, stats::sd) > 0])

  # Same as above but for wy.mat:
  wy.mat <- as.matrix(wy.mat[, apply(wy.mat, 2, stats::sd) > 0])

  # Same as above but for wXs.mat:
  wXs.mat <- as.matrix(wXs.mat[, apply(wXs.mat, 2, stats::sd) > 0])

  # Same as above but for wys.mat:
  wys.mat <- as.matrix(wys.mat[, apply(wys.mat, 2, stats::sd) > 0])

  # Gives the colnames (i.e., covariate names) from `X` (`y`) to rows of
  #   `wX.mat` (`wy.mat`) and `wXs.mat` (`wys.mat`).
  rownames(wXs.mat) <- rownames(wX.mat) <- colnames(X)
  rownames(wys.mat) <- rownames(wy.mat) <- colnames(y)

  # Similar to above two lines of code, but all the `Z` matrices.
  rownames(Z.mat) <- rownames(ZX.mat) <- rownames(Zy.mat) <- rownames(X)

  # If there are shared covariates, regresses `Z.mat` on `X.c` and assigns
  #   the resulting `lm` object to `beta.out`.
  if (length(X.c) > 0) {
    beta.out <- stats::lm(Z.mat ~ X.c)
  } else {
    beta.out <- NULL
  }

  output <- list(
    "z" = Z.mat,
    "z.X" = ZX.mat,
    "z.y" = Zy.mat,
    "w.Xs" = wXs.mat,
    "w.ys" = wys.mat,
    "w.X" = wX.mat,
    "w.y" = wy.mat,
    "beta.z" = beta.out,
    "lz" = lz,
    "lz.X" = lz.X,
    "lz.y" = lz.y,
    "bic" = bic.sort,
    "proportionX" = proportionX
  )

  return(output)
}


