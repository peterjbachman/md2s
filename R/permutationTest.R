
## ---------------------------------------------------------------------------------------------------
##
## 		Replication Files: Create Data Figure 3 (Panel A)
##		"Scaling Data from Multiple Sources"
##		Authors: Ted Enamorado, Gabriel Lopez-Moctezuma, Marc Ratkovic
##
## ---------------------------------------------------------------------------------------------------
source('./R/00_MD2S.R')
if(exists("out.all")){rm(out.all)}

kX.num <- 100
n <- 50
ky <- 40
## the line below comes from Run.R for the simulation
## Number of simulations; permutations; and bootstap samples
nsims <- nperm <- nboot <- 10 # originally 200, too long to run :)


##  tmp code setup
# kX = 100
# j = 1

##Generate Data
results.all<-NULL
for(j in 1:nsims) {
  for(kX in c(kX.num)){

    ## -------------------------
    ## Create z.s y z2.s
    ## -------------------------
    z1 <- rnorm(n); z2 <- rnorm(n)
    zs.true <- cbind(z1, z2)
    zX.true <- cbind(rnorm(n), rnorm(n), rnorm(n))
    zy.true <- cbind(rnorm(n), rnorm(n))

    bsX.true <- cbind(rnorm(kX), rnorm(kX))
    bX.true <- cbind(rnorm(kX), rnorm(kX), rnorm(kX))
    bsy.true <- cbind(rnorm(ky), rnorm(ky))
    by.true <- cbind(rnorm(ky), rnorm(ky))

    dw.z <- diag(c(2, 1))
    dw.y <- diag(c(4, 2))
    dw.X <- diag(c(4, 2, 2))

    X1 <- zs.true %*% dw.z %*% t(bsX.true) +  zX.true %*% dw.X %*% t(bX.true) + 2 * matrix(rnorm(n * kX), nr = n)
    y1 <- zs.true %*% dw.z %*% t(bsy.true) +  zy.true %*% dw.y %*% t(by.true) + 2 * matrix(rnorm(n * ky), nr = n)


    # Notice here b1 is a MD2S obejct #
    b1 <- MD2S(X=(X1),y=(y1), dim = 5)

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

    stopCluster(cl)

    lz.run <- lz.X.run <- lz.y.run <- list()
    for(i in 1:length(results)) {
      lz.run[[i]] <- results[[i]]$lz
      lz.X.run[[i]] <- results[[i]]$lz.X
      lz.y.run[[i]] <- results[[i]]$lz.y
    }

    lz.run2  <- do.call('rbind', lz.run)
    lz.X.run2 <- do.call('rbind', lz.X.run)
    lz.y.run2 <- do.call('rbind', lz.y.run)

    pz1<-pz2<-pz2.y<-pz2.X<-NULL
    dims <- 5

    for(i in 1:dims){
      pz2[i]<- mean(b1$lz[i] < lz.run2[,i])
      pz2.y[i] <- mean(b1$lz.y[i] < lz.y.run2[,i])
      pz2.X[i] <- mean(b1$lz.X[i] < lz.X.run2[,i])
    }

    pvals.s <- pz2
    pvals.X <- pz2.X
    pvals.y <- pz2.y

    pvals.s
    pvals.X
    pvals.y

    scenario2 <- 1
    results.curr<-c(n,kX,scenario2,pvals.s,pvals.X,pvals.y)

    names(results.curr)<-c("N","kX","scen2",paste("ps_",1:5,sep=""),paste("pX_",1:5,sep=""),paste("py_",1:5,sep=""))

    results.all <- results.curr

    ran.name <- paste0("./results_Panel_A/output_",round(runif(1),10)*1e10, "_Panel_A.RData")
    save(results.all, file = ran.name)
  }
}
