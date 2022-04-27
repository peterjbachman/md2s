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

# SECTION 1

MD2S <- function(X,y,X.s=NULL,X.X=NULL, X.y=NULL,init="svd",dim=1,sim=FALSE,tol=1e-6,BIC=FALSE){
  if(BIC==FALSE) bic.sort<-NULL
  X.c<-X.s
  n<-nrow(X)
  ZX.mat<-Zy.mat<-Z.mat<-as.matrix(rep(1,n))

  wXs.mat<-wX.mat<-rep(1,ncol(X))
  wys.mat<-wy.mat<-rep(1,ncol(y))
  ZX.mat<-Zy.mat<-Z.mat<-as.matrix(rep(1,n))
  proportionX<-lz<-lz.X<-lz.y<-NULL
  X1<-X-make.int(X)
  y1<-y-make.int(y)
  for(i in 1:dim){
    cat("----  Fitting Dimension ", i, "----  \n")
    X1<-fastres(X1,cbind(Z.mat,ZX.mat));y1<-fastres(y1,cbind(Z.mat,Zy.mat))
    X1<-t(fastres(t(X1),cbind(wX.mat,wXs.mat)))
    y1<-t(fastres(t(y1),cbind(wy.mat,wys.mat)))
    X1<-X1-make.int(X1)
    y1<-y1-make.int(y1)
    b0<-MD2S_inner(X0=X1,y0=y1,X.c.0=X.c,X.X.0=X.X, X.y.0=X.y,tol0=tol)#,burnin=0,gibbs=100,thin=1)
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

MD2S_inner <- function(X0,y0,X.c.0=NULL,X.X.0=NULL, X.y.0=NULL,init0="svd",sim0=FALSE,tol0=tol){

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

## Creates integers, but based on what?
make.int <- function(X) {
  int1 <- rep(1, nrow(X)) %*% t(colMeans(X))
  int2 <- rowMeans(X) %*% t(rep(1, ncol(X)))
  int1 + int2 - mean(int1 + int2) + mean(X)
}

## Create residuals faster?
## TODO Add library MASS to package
fastres <- function(x, z) {
  z <- cbind(1, z)
  # If there are fewer rows than columns
  if (nrow(z) <= ncol(z)) fits <- z %*% MASS::ginv(t(z) %*% z) %*% (t(z) %*% x)
  # If there are more columns than rows
  if (nrow(z) > ncol(z)) fits <- z %*% (MASS::ginv(t(z) %*% z) %*% t(z) %*% x)

  # Create residuals
  res <- x - fits
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
