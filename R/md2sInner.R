#I think that this is getting some inner product space. It's using methods to compute solutions to an underdetermined system.
md2sInner <- function(X0, # Must be the double-centered/scaled matrix derived from X.
                       y0, # Must be the double-centered/scaled matrix derived from y.
                       X.c.0 = NULL, # covariates associated with shared subspace.
                       X.X.0 = NULL, # covariates associated with X subspace.
                       X.y.0 = NULL, # covariates associated with y subspace.
                       init0 = "svd", # singular value decomposition (no other options?)
                       sim0 = FALSE, # Indicates whether we're using simulated data (`TRUE` not an option?)
                       tol0 = tol # Convergence/iteration tolerance parameter.
) {
  X <- X0
  y <- y0
  X.c <- X.c.0
  init <- init0
  sim <- sim0

  tol <- tol0

  X.X <- X.X.0
  X.y <- X.y.0

  X.c <- cleanup(X.c)
  X.X <- cleanup(X.X)
  X.y <- cleanup(X.y)


  ## Declare and initialize, same as above.
  n <- nrow(X)
  X1 <- X
  y1 <- y
  z.x <- z.y <- z <- rep(1, n)

  #Okay, so I THINK this is getting the generalized inverse
  XXprime <- X %*% t(X)
  yyprime <- y %*% t(y)
  #If dim of these is non-0--so I think if the system is overdetermined, then find generallized inverse.
  if (length(X.c) > 0) hat.Xc <- MASS::ginv(t(X.c) %*% X.c) %*% t(X.c)
  if (length(X.X) > 0) hat.XX <- MASS::ginv(t(X.X) %*% X.X) %*% t(X.X)
  if (length(X.y) > 0) hat.Xy <- MASS::ginv(t(X.y) %*% X.y) %*% t(X.y)
  #So this is telling us which way to solve the system.
  if (init == "svd") {
    ##Whenever we have FALSE passing through, we run singular value decomp to solve. This, I THINK, gets the left singular value decomp.
    if (FALSE) {
      z.x <- svd(X1, nu = 1)$u[, 1]
      z.y <- svd(y1, nu = 1)$u[, 1]
    }
    ##Okay, yes. Whenever X1 has more cols then rows--underdetermined--make z.x solution from psudeoinverse. Otherwise, use my.norm over the
    ##product of X1 and the pseudoinverse of the crossproduct of X1. I think that this comes from
    if (nrow(X1) < ncol(X1)) z.x <- irlba::irlba(X1, nu = 1, nv = 1)$u[, 1] else z.x <- my.norm(X1 %*% irlba::irlba(crossprod(X1), nu = 1, nv = 1)$v[, 1])
    ##Same as above but for y1.
    if (nrow(y1) < ncol(y1)) z.y <- irlba::irlba(y1, nu = 1, nv = 1)$u[, 1] else z.y <- my.norm(y1 %*% irlba::irlba(crossprod(y1), nu = 1, nv = 1)$v[, 1])

    z <- (z.x + z.y) / 2
    z <- my.norm(z)
    z.x <- my.norm(z.x)
    z.y <- my.norm(z.y)
  }

  loglik <- 0

  # SECTION 2
  ##I think this measures the correlation between true and estimated values over 1000 simulations per combo of N and K_1. I'm confident that this is making the graphics that we see in the appendix.
  for (i in 1:1000) {

    X1 <- X - (z.x %*% t(z.x)) %*% X
    y1 <- y - (z.y %*% t(z.y)) %*% y

    eig.form <- t(X1) %*% y1
    eig.form2 <- t(y1) %*% X1

    if (FALSE) {
      svd1 <- svd(eig.form, nu = 1)
      svd2 <- svd(eig.form2, nu = 1)
    }
    svd1 <- irlba::irlba(eig.form, nu = 1, nv = 1)
    svd2 <- irlba::irlba(eig.form2, nu = 1, nv = 1)
    w1 <- (svd1$u[, 1])
    w2 <- (svd2$u[, 1])
    z.freq1 <- my.norm(X1 %*% w1)
    z.freq2 <- my.norm(y1 %*% w2)
    z.freq2 <- z.freq2 * sign(stats::cor(z.freq1, z.freq2))

    ## Initialize least squares estimates
    if (length(X.c) > 0) {
      z.freq3 <- as.vector(X.c %*% (hat.Xc %*% z))
      z.freq3 <- z.freq3 * stats::cor(z.freq2, z.freq3)
      z.freq3 <- my.norm(z.freq3)
    }

    if (length(X.X) > 0) {
      z.freq3X <- as.vector(X.X %*% (hat.XX %*% z.x))
      z.freq3X <- z.freq3X * stats::cor(z.freq1, z.freq3X)
      z.freq3X <- my.norm(z.freq3X)
    }

    if (length(X.y) > 0) {
      z.freq3y <- as.vector(X.y %*% (hat.Xy %*% z.y))
      z.freq3y <- z.freq3y * stats::cor(z.freq2, z.freq3y)
      z.freq3y <- my.norm(z.freq3y)
    }

    z <- my.norm(z)

    cor.last <- check.cor(z)

    z.last <- z

    # Can be streamlined if we use the
    alpha.min <- stats::optimize(alpha.func1, lower = -5, upper = 5, maximum = TRUE)$max
    p1 <- exp(alpha.min)
    p1 <- p1 / (1 + p1)
    z <- scale(p1 * z.freq1 + (1 - p1) * z.freq2)
    z <- my.norm(z)
    p1.f <- p1 / (1 + p1)

    ## Update with covariates
    if (length(X.c) > 0) {
      lm.z <- stats::lm(z ~ z.freq3)
      z.fit3 <- my.norm(lm.z$fit)
      z.res3 <- my.norm(lm.z$res)
      alpha.min <- stats::optimize(alpha.func2, lower = -5, upper = 5, maximum = TRUE)$max
      p1 <- exp(alpha.min)
      p1 <- p1 / (1 + p1)
      z <- scale(p1 * z.fit3 + (1 - p1) * z.res3)
      z <- my.norm(z)
    }

    if (length(X.X) > 0) {
      lm.z <- stats::lm(z.x ~ z.freq3X)
      z.fit3 <- my.norm(lm.z$fit)
      z.res3 <- my.norm(lm.z$res)
      alpha.min <- stats::optimize(alpha.func2.X, lower = -5, upper = 5, maximum = TRUE)$max
      p1 <- exp(alpha.min)
      p1 <- p1 / (1 + p1)
      z.x <- scale(p1 * z.fit3 + (1 - p1) * z.res3)
      z.x <- my.norm(z.x)
    }

    if (length(X.y) > 0) {
      lm.z <- stats::lm(z.y ~ z.freq3y)
      z.fit3 <- my.norm(lm.z$fit)
      z.res3 <- my.norm(lm.z$res)
      alpha.min <- stats::optimize(alpha.func2.y, lower = -5, upper = 5, maximum = TRUE)$max
      p1 <- exp(alpha.min)
      p1 <- p1 / (1 + p1)
      z.y <- scale(p1 * z.fit3 + (1 - p1) * z.res3)
      z.y <- my.norm(z.y)
    }


    loglik.last <- loglik
    lz <- sum((t(X1) %*% z)^2) + sum((t(y1) %*% z)^2)
    lz.x <- sum((t(X1) %*% z.x)^2)
    lz.y <- sum((t(y1) %*% z.y)^2)
    loglik <- lz + lz.x + lz.y

    if (nrow(X1) < ncol(X1)) {
      loglik <- t(z) %*% ((X1 %*% t(X1)) %*% (y1 %*% t(y1))) %*% z
    } else {
      loglik <- (t(z) %*% X1) %*% t(X1) %*% y1 %*% (t(y1) %*% z)
    }
    loglik <- as.vector(loglik)

    if (i %% 50 == 0) {
      cat("  Iteration ", i, "\n")
      cat("  Current log-likelihood bound:", round(loglik, 4), "\n")
    }
    if (i > 1) {
      if (abs(loglik - loglik.last) / loglik.last < tol | (loglik < loglik.last & i > 10)) {
        cat("## Convergence after ", i, "iterations ## \n\n")
        break
      }
    }

    Xlessz <- fastres(X, z)
    ylessz <- fastres(y, z)

    z.x.try <- my.norm(irlba::irlba(Xlessz, nu = 1, nv = 1)$u[, 1])
    z.y.try <- my.norm(irlba::irlba(ylessz, nu = 1, nv = 1)$u[, 1])

    if (i == 1) {
      z.x <- z.x.try
      z.y <- z.y.try
    }
    if (i > 1) {
      z.x <- z.x.try * sign(stats::cor(z.x, z.x.try))
      z.y <- z.y.try * sign(stats::cor(z.y, z.y.try))
    }
    z.x <- my.norm(z.x)
    z.y <- my.norm(z.y)
    z <- my.norm(z)
    w.Xs <- my.norm(t(z) %*% X1)
    w.ys <- my.norm(t(z) %*% y1)

    if (sim) print(c(stats::cor(z, z.true), stats::cor(z.x, zX.true), stats::cor(z.y, zy.true)))
  }


  Xc.out <- NULL
  if (length(X.c) > 0) {
    Xc.out <- stats::lm(z ~ X.c)$coef[-1]
  }

  output <- list("z" = z, "z.X" = z.x, "z.y" = z.y, "w.Xs" = w.Xs, "w.ys" = w.ys, "beta.z" = Xc.out, "proportionX" = p1.f)

  return(output)
}
