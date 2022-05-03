## ---------------------------------------------------------------------------------------------------
##
## 		Replication Files: MD2S function
##		"Scaling Data from Multiple Sources"
##		Authors: Ted Enamorado, Gabriel Lopez-Moctezuma, Marc Ratkovic
##
## ---------------------------------------------------------------------------------------------------

## Required libraries:
library(MASS)
library(irlba)

# SECTION 1 - wowsa this is gonna be fun! ::UpsideDownSmileyFaceEmoji::

MD2S <- function( # List of arguments and descriptions below:
  X, # N x K_1 dataset m = 1 of realized outcomes.
    # DESCRIPTION: first dataset of realized outcomes for `N` observations and
    #   `K_1` covariates. Note that `K_1` & `K_2` can differ.
    # EXAMPLE: roll-call votes for `N` legislators with `K_1` covariates.
  
  y, # N x K_2 dataset m = 2 of realized outcomes.
    # DESCRIPTION: first dataset of realized outcomes for `N` observations and
    #   `K_2` covariates. Note that `K_1` & `K_2` can differ.
    # EXAMPLE: floor speech text for `N` legislators with `K_2` covariates.
  
  X.s = NULL, # Optional covariates for estimating scaled locations in shared subspace.
    # DESCRIPTION: dataset of covariates that explain observations' scaled 
    #   locations in shared subspace.
    # EXAMPLE: covariates that explain legislators' general ideological placements.
  
  X.X = NULL, # Optional covariates for estimating scaled locations in idiosyncratic subspace of X.
    # DESCRIPTION: dataset of covariates that explain observations' scaled 
    #   locations in idiosyncratic subspace of X.
    # EXAMPLE: covariates that explain legislators' ideological placements as
    #   determined by roll call votes.
  
  X.y = NULL, # Optional covariates for estimating scaled locations in idiosyncratic subspace of y.
    # DESCRIPTION: dataset of covariates that explain observations' scaled 
    #   locations in idiosyncratic subspace of y.
    # EXAMPLE: covariates that explain legislators' ideological placements as
    #   determined by floor speech text.
  
  init = "svd", # Initialize??? 
    # NOTE: Default is SVD = Singular Value Decomposition of a matrix.
    # NOTE: Function permits no other argument besides default "svd".
    # [VERIFICATION CHECK NEEDED]
  
  dim = 1, # Number of fitted dimensions (from most- to least-explanatory).
    # (???) Not entirely sure if this is correct interpretation, check (???)
    # [VERIFICATION CHECK NEEDED]
  
  sim = FALSE, # `TRUE` if using simulated data (???)
    # (???) May need to have "true" values for z.true, zX.true, zy.true in environment. (???)
    # Perhaps this part should be removed given that z.true, zX.true, & zy.true 
    #   are not created in the function, so users cant run `sim = TRUE` unless
    #   they have objects named z.true, zX.true, zy.true in their local environment.
    #   Alternatively, we could force users to provide true values if `sim = TRUE`.
    # [VERIFICATION CHECK NEEDED]
  
  tol = 1e-6, # TOL = convergence tolerance parameter. 
      # DESCRIPTION: criteria for iteration termination. Iterations terminate
      #   if difference between current & prior iteration < `TOL` value.
    # [VERIFICATION CHECK NEEDED] (though I'm pretty sure this is right)
  
  BIC = FALSE # (??BIC = Bayesian Information Criterion??) 
    # Very unsure on this. This argument `BIC` apparently doesn't do anything
    #   in the function besides appearing as an element in the `output` list.
    # [VERIFICATION CHECK NEEDED]
  ) {
    
    # Sets `bic.sort` to `NULL` if the `BIC` argument equals `FALSE`, which sets the "bic" element of the `output` list to `NULL`
    if(BIC==FALSE) bic.sort <- NULL
    
    # Renames optional covariates for explaining scaled locations in the shared subspace.
    X.c <- X.s 
    
    # Sets `n` as the number of observations in X, the first dataset (m = 1) of realized outcomes.
    n <- nrow(X)
    
    # Creates three `Z` column vectors of 1s with length = number of rows in X.
    #   (???) Correspond to the shared (Z_S) & idiosyncratic (Z_{(M)}) subspaces.
    #     (from paper): Z_S contains latent locations in the shared subspace in columns for each `dim` dimensions (in this case, 1).
    #     (from paper): Z_{(M)} contains latent locations in the idiosyncratic subspace for `dim` latent dimensions (in this case, 1).
    ZX.mat <- Zy.mat <- Z.mat <- as.matrix(rep(1, n))

    
    # Creates two (row?) vectors of 1s with length = number of columns in X.
    #   (from paper): W_{(M)} is matrix of shared factors for the shares subspace for the dataset Y_{(M)}.
    wXs.mat <- wX.mat <- rep(1,ncol(X))
    
    # Creates two (row?) vectors of 1s with length = number of columns in y.
    #   (from paper): W_{(M)} is matrix of shared factors for the shares ubspace for the dataset Y_{(M)}.
    wys.mat <- wy.mat <- rep(1, ncol(y))
    
    # REPEAT of third to last line of code.
    #   Delete unnecessary repeat?
    ZX.mat <- Zy.mat <- Z.mat <- as.matrix(rep(1, n))
    
    # Creates `NULL` objects (??for later use??).
    #   Potential purpose: proportionX may be intended to measure how much larger the `X` dataset is than the `y` dataset. It does nothing here, though.
    proportionX <- lz <- lz.X <- lz.y <- NULL
    
    # Double-center `X1` and `y1`:
    #   "we preprocess the matrices by double-centering them, so that the row-mean, 
    #   column-mean, and grand mean is zero" (pg. 216 of paper)
    # BACKGROUND: the `make.int` function is defined far below. 
    # `make.int` standardizes the rows (columns) of a matrix such that each
    #   row (column) of the matrix has the same standard deviation & mean.
    X1 <- X - make.int(X)
    y1 <- y - make.int(y)
    
    #
    for(i in 1:dim){
      
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
      y1<-t(fastres(t(y1),cbind(wy.mat,wys.mat)))
      
      # double-centers new X1:
      #   Again, means of new result is basically closer to zero.
      X1<-X1 - make.int(X1)
      
      # Same as directly above but for new y1:
      y1<-y1-make.int(y1)
      
      # (???) WHAT DOES `MD2S_inner` FUNCTION DO (???)
      b0 <- MD2S_inner(
        X0 = X1, # Double-centered/scaled, derived from X.
        y0 = y1, # Double-centered/scaled, derived from y1.
        X.c.0 = X.c, # covariates associated with shared subspace.
        X.X.0 = X.X, # Covariates associated with X subspace.
        X.y.0 = X.y, # covariates associated with y subspace.
        tol0 = tol   # tolerance parameter.
      )
      #,burnin=0,gibbs=100,thin=1) - ORIGINAL COMMENT! Potentially additional arguments they planned to provide??
      
      Z.mat<-cbind(Z.mat,b0$z)
      
      Zy.mat<-cbind(Zy.mat,b0$z.y)
      
      ZX.mat<-cbind(ZX.mat,b0$z.X)
      
      wXs.mat<-cbind(wXs.mat,my.norm(t(b0$z%*%X1)))
    wys.mat<-cbind(wys.mat,my.norm(t(b0$z%*%y1)))
    wX.mat<-cbind(wX.mat,my.norm(t(b0$z.X%*%X1)))
    wy.mat<-cbind(wy.mat,my.norm(t(b0$z.y%*%y1)))
    proportionX[i]<-b0$pr

    if(nrow(X1)<max(ncol(X1),ncol(y1))){
      lz[i]<-t(b0$z)%*%(X1%*%t(X1))%*%(y1%*%t(y1))%*%b0$z
      } else {
        lz[i]<-((t(b0$z)%*%X1)%*%t(X1))%*%(y1%*%(t(y1)%*%b0$z))
      }
      lz.X[i]<-sum((t(X1)%*%b0$z.X)^2)
      lz.y[i]<-sum((t(y1)%*%b0$z.y)^2)
    }

    Z.mat<-as.matrix(Z.mat[,apply(Z.mat,2,sd)>0])
    Zy.mat<-as.matrix(Zy.mat[,apply(Zy.mat,2,sd)>0])
    ZX.mat<-as.matrix(ZX.mat[,apply(ZX.mat,2,sd)>0])
    wX.mat<-as.matrix(wX.mat[,apply(wX.mat,2,sd)>0])
    wy.mat<-as.matrix(wy.mat[,apply(wy.mat,2,sd)>0])
    wXs.mat<-as.matrix(wXs.mat[,apply(wXs.mat,2,sd)>0])
    wys.mat<-as.matrix(wys.mat[,apply(wys.mat,2,sd)>0])

    rownames(wXs.mat)<-rownames(wX.mat)<-colnames(X)
    rownames(wys.mat)<-rownames(wy.mat)<-colnames(y)

    rownames(Z.mat)<-rownames(ZX.mat)<-rownames(Zy.mat)<-rownames(X)

    if(length(X.c)>0) beta.out<-lm(Z.mat~X.c) else beta.out<-NULL
    output<-list("z"=Z.mat,"z.X"=ZX.mat,"z.y"=Zy.mat,"w.Xs"=wXs.mat,"w.ys"=wys.mat,
               "w.X"=wX.mat,"w.y"=wy.mat, "beta.z"=beta.out,"lz"=lz,"lz.X"=lz.X,"lz.y"=lz.y,"bic"=bic.sort,"proportionX"=proportionX)

    return(output)
}

MD2S_inner <- function(
  X0, # Must be the double-centered/scaled matrix derived from X.
  y0, # Must be the double-centered/scaled matrix derived from y.
  X.c.0 = NULL, # covariates associated with shared subspace.
  X.X.0 = NULL, # covariates associated with X subspace.
  X.y.0 = NULL, # covariates associated with y subspace.
  init0 = "svd", # singular value decomposition (no other options?)
  sim0 = FALSE, # Indicates whether we're using simulated data (`TRUE` not an option?)
  tol0 = tol # Convergence/iteration tolerance parameter.
  ) {

    X<-X0;y<-y0;X.c<-X.c.0;init<-init0;sim<-sim0;

    tol<-tol0

    X.X<-X.X.0;X.y<-X.y.0

    my.norm<-function(x) {x<-as.vector(x);x<-x-mean(x);(x/sum(x^2)^.5)}
    cleanup<-function(X.c){
      if(length(X.c)>0){
        X.c<-apply(X.c,2,FUN=function(x) x-mean(x))
        X.c<-X.c[,colMeans(X.c^2)>1e-10]
      }
      X.c
    }
    X.c<-cleanup(X.c)
    X.X<-cleanup(X.X)
    X.y<-cleanup(X.y)


    ##Declare and initialize
    n<-nrow(X)
    X1<-X;y1<-y
    z.x<-z.y<-z<-rep(1,n)

    XXprime<-X%*%t(X)
    yyprime<-y%*%t(y)

    if(length(X.c)>0) hat.Xc<-ginv(t(X.c)%*%X.c)%*%t(X.c)
    if(length(X.X)>0) hat.XX<-ginv(t(X.X)%*%X.X)%*%t(X.X)
    if(length(X.y)>0) hat.Xy<-ginv(t(X.y)%*%X.y)%*%t(X.y)
  
    if(init=="svd"){
      if(FALSE){
        z.x<-svd(X1,nu=1)$u[,1]
        z.y<-svd(y1,nu=1)$u[,1]
      }
      if(nrow(X1)<ncol(X1)) z.x<-irlba(X1,nu=1,nv=1)$u[,1] else z.x<-my.norm(X1%*%irlba(crossprod(X1),nu=1,nv=1)$v[,1])
      if(nrow(y1)<ncol(y1)) z.y<-irlba(y1,nu=1,nv=1)$u[,1] else z.y<-my.norm(y1%*%irlba(crossprod(y1),nu=1,nv=1)$v[,1])
  
      z<-(z.x+z.y)/2
      z<-my.norm(z);z.x<-my.norm(z.x);z.y<-my.norm(z.y)
  
    }
  
    loglik<-0
  
  # SECTION 2
  
    for(i in 1:1000){
  
      rm.zX<-function(x) fastres(x,z.x)
      rm.zy<-function(x) fastres(x,z.y)
      rm.z<-function(x) fastres(x,z)
  
      X1<-X-(z.x%*%t(z.x))%*%X
      y1<-y-(z.y%*%t(z.y))%*%y
  
      eig.form<-t(X1)%*%y1
      eig.form2<-t(y1)%*%X1
  
      if(FALSE){
        svd1<-svd(eig.form,nu=1)
        svd2<-svd(eig.form2,nu=1)
      }
      svd1<-irlba(eig.form,nu=1,nv=1)
      svd2<-irlba(eig.form2,nu=1,nv=1)
      w1<-(svd1$u[,1])
      w2<-(svd2$u[,1])
      z.freq1<-my.norm(X1%*%w1)
      z.freq2<-my.norm(y1%*%w2)
      z.freq2<-z.freq2*sign(cor(z.freq1,z.freq2))
    
      ##Initialize least squares estimates
      if(length(X.c)>0){
        z.freq3<-as.vector(X.c%*%(hat.Xc%*%z))
        z.freq3<-z.freq3*cor(z.freq2,z.freq3)
        z.freq3<-my.norm(z.freq3)
      }
  
      if(length(X.X)>0){
        z.freq3X<-as.vector(X.X%*%(hat.XX%*%z.x))
        z.freq3X<-z.freq3X*cor(z.freq1,z.freq3X)
        z.freq3X<-my.norm(z.freq3X)
      }
  
      if(length(X.y)>0){
        z.freq3y<-as.vector(X.y%*%(hat.Xy%*%z.y))
        z.freq3y<-z.freq3y*cor(z.freq2,z.freq3y)
        z.freq3y<-my.norm(z.freq3y)
      }
  
      z<-my.norm(z)
  
      alpha.func<-function(x,z1,z2){
        p1<-exp(x)/(1+exp(x))
        check.cor(p1*z1+(1-p1)*z2)
      }
  
      alpha.func.X<-function(x,z1,z2){
        p1<-exp(x)/(1+exp(x))
        check.cor.X(p1*z1+(1-p1)*z2)
      }
  
      alpha.func.y<-function(x,z1,z2){
        p1<-exp(x)/(1+exp(x))
        check.cor.y(p1*z1+(1-p1)*z2)
      }
  
      alpha.func1<-function(x) alpha.func(x,z.freq1,z.freq2)
      alpha.func2<-function(x) alpha.func(x,z.fit3,z.res3)
      alpha.func2.X<-function(x) alpha.func.X(x,z.fit3,z.res3)
      alpha.func2.y<-function(x) alpha.func.y(x,z.fit3,z.res3)
  
      check.cor<-function(z.run){
        z.run<-my.norm(z.run-mean(z.run))
        if(nrow(X)<ncol(X)) wX<-as.vector(XXprime%*%z.run) else wX<-as.vector(X1%*%(t(X1)%*%z.run))
        if(nrow(y)<ncol(y)) wy<-as.vector(yyprime%*%z.run) else wy<-as.vector(y1%*%(t(y1)%*%z.run))
        cov(wX,wy)
      }
  
      check.cor.X<-function(z.run){
        z.run<-my.norm(z.run-mean(z.run))
        if(nrow(X)<ncol(X)) wX<-as.vector(XXprime%*%z.run) else wX<-as.vector(X%*%(t(X)%*%z.run))
        var(wX)
      }
  
      check.cor.y<-function(z.run){
        z.run<-my.norm(z.run-mean(z.run))
        if(nrow(y)<ncol(y)) wy<-as.vector(yyprime%*%z.run) else wy<-as.vector(y%*%(t(y)%*%z.run))
        var(wy)
      }
      cor.last<-check.cor(z)
  
      z.last<-z
      alpha.min<-optimize(alpha.func1,lower=-5,upper=5,maximum=TRUE)$max
      p1<-exp(alpha.min)
      p1<-p1/(1+p1)
      z<-scale(p1*z.freq1+(1-p1)*z.freq2)
      z<-my.norm(z)
      p1.f <-p1/(1+p1)
  
      ##Update with covariates
      if(length(X.c)>0){
        lm.z<-lm(z~z.freq3)
        z.fit3<-my.norm(lm.z$fit)
        z.res3<-my.norm(lm.z$res)
        alpha.min<-optimize(alpha.func2,lower=-5,upper=5,maximum=TRUE)$max
        p1<-exp(alpha.min)
        p1<-p1/(1+p1)
        z<-scale(p1*z.fit3+(1-p1)*z.res3)
        z<-my.norm(z)
      }
  
      if(length(X.X)>0){
        lm.z<-lm(z.x~z.freq3X)
        z.fit3<-my.norm(lm.z$fit)
        z.res3<-my.norm(lm.z$res)
        alpha.min<-optimize(alpha.func2.X,lower=-5,upper=5,maximum=TRUE)$max
        p1<-exp(alpha.min)
        p1<-p1/(1+p1)
        z.x<-scale(p1*z.fit3+(1-p1)*z.res3)
        z.x<-my.norm(z.x)
      }
  
      if(length(X.y)>0){
        lm.z<-lm(z.y~z.freq3y)
        z.fit3<-my.norm(lm.z$fit)
        z.res3<-my.norm(lm.z$res)
        alpha.min<-optimize(alpha.func2.y,lower=-5,upper=5,maximum=TRUE)$max
        p1<-exp(alpha.min)
        p1<-p1/(1+p1)
        z.y<-scale(p1*z.fit3+(1-p1)*z.res3)
        z.y<-my.norm(z.y)
      }
  
  
      loglik.last<-loglik
      lz<-sum((t(X1)%*%z)^2)+sum((t(y1)%*%z)^2)
      lz.x<-sum((t(X1)%*%z.x)^2)
      lz.y<-sum((t(y1)%*%z.y)^2)
      loglik<-lz+lz.x+lz.y

     if(nrow(X1)<ncol(X1)) {
        loglik<-t(z)%*%((X1%*%t(X1))%*%(y1%*%t(y1)))%*%z
      } else {
        loglik<-(t(z)%*%X1)%*%t(X1)%*%y1%*%(t(y1)%*%z)
      }
      loglik<-as.vector(loglik)
  
      if(i%%50==0)	{
        cat("  Iteration ", i, "\n")
        cat("  Current log-likelihood bound:", round(loglik,4), "\n")
  
  
      }
      if(i>1)  {
        if(abs(loglik-loglik.last)/loglik.last < tol| (loglik<loglik.last & i>10) ){
          cat("## Convergence after ", i, "iterations ## \n\n")
          break
        }
      }

      Xlessz<-rm.z(X)
      ylessz<-rm.z(y)
  
      z.x.try<-my.norm(irlba(Xlessz,nu=1,nv=1)$u[,1])
      z.y.try<-my.norm(irlba(ylessz,nu=1,nv=1)$u[,1])
  
      if(i==1) {z.x<-z.x.try; z.y<-z.y.try}
      if(i>1) {
        z.x<-z.x.try*sign(cor(z.x,z.x.try))
        z.y<-z.y.try*sign(cor(z.y,z.y.try))
      }
      z.x<-my.norm(z.x);z.y<-my.norm(z.y);z<-my.norm(z)
      w.Xs<-my.norm(t(z)%*%X1);w.ys<-my.norm(t(z)%*%y1)

      if(sim) print(c(cor(z,z.true),cor(z.x,zX.true),cor(z.y,zy.true)))

    }


    Xc.out<-NULL
    if(length(X.c)>0){Xc.out<-lm(z~X.c)$coef[-1]}

    output<-list("z"=z,"z.X"=z.x,"z.y"=z.y,"w.Xs"=w.Xs,"w.ys"=w.ys,"beta.z"=Xc.out,"proportionX"=p1.f)

    return(output)
}

# SECTION 3

## Checks the correlation and outputs it?
## TODO Add inputs for x1 and y1 for my sanity
check.cor <- function(z.run) {
  wX <- as.vector(((X1) %*% t(X1)) %*% z.run)
  wy <- as.vector(((y1) %*% t(y1)) %*% z.run)
  cov(wX, wy)
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
## TODO Add library MASS to package
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
my.norm <- function(x) {
  x <- as.vector(x)
  x <- x - mean(x)
  (x / sum(x^2)^.5)
}

## No clue yet
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
  qnorm(edf)
}
