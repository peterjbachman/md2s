## ---------------------------------------------------------------------------------------------------
##
## 		Replication Files: Create Data Figure 3 (Panel A)
##		"Scaling Data from Multiple Sources"
##		Authors: Ted Enamorado, Gabriel Lopez-Moctezuma, Marc Ratkovic
## 
## ---------------------------------------------------------------------------------------------------
source('./00_MD2S.R')
if(exists("out.all")){rm(out.all)}

## ---------------------------------------------------------------------------------------------------
## Set up see page 224 in the paper 
## Panel A: N = 50, K(2) = 100
## Panel B: N = 100, K(2) = 100
## Panel C: N = 100, K(2) = 1000
## For each of the above-mentioned settings, 1,000 total simulations with
## 1,000 permuted datasets per simulation were used to estimate the p value.
## 
## The numbers should be made available for users to define in the function parameters. 
## Set default to the ones below: 
## ---------------------------------------------------------------------------------------------------
kX.num <- 100 # k(2) = 100 for Panel A 
n <- 50 # N = 50 for Panel A 
ky <- 40 # k(1) is held at 40
# k(1) and k(2) are columns in the matrices 


## ---------------------------------------------------------------------------------------------------
## the line below comes from Run.R for the simulation
## Number of simulations; permutations; and bootstap samples
## ---------------------------------------------------------------------------------------------------
nsims <- nperm <- nboot <- 200 # originally 200, too long to run :) 



## ---------------------------------------------------------------------------------------------------
## Generate Data
## ---------------------------------------------------------------------------------------------------
results.all <- NULL # create list to store outputs 

# outter loop runs nsims <- 10 times, originally 200 
for(j in 1:nsims) {
  # inner loop basically runs once with kX.num <- 100 set above 
  # where kX = kX.num
  for(kX in c(kX.num)){
    
    ## ---------------------------------------------------------------------------
    ## Create z.s y z2.s
    ## See page 221: 
    ## All systematic factors Z, W, B are drawn from a standard normal 
    ## Y(1) has 2 idiosyncratic dimensions --> Y
    ## Y(2) has 3 dimensions --> X
    ## ---------------------------------------------------------------------------
    
    # z1 <- rnorm(n) # 50 numbers generated from normal distribution with mean = 0 sd = 1
    # z2 <- rnorm(n) # 50 numbers generated from normal distribution with mean = 0 sd = 1
    # zs.true <- cbind(z1, z2) # cbind 2 cols z1, z2 # dim 50, 2
    
    zs.true <- cbind(rnorm(n), rnorm(n)) # I combined the previous 3 lines to be more consistent
    zX.true <- cbind(rnorm(n), rnorm(n), rnorm(n)) # 50, 3
    zy.true <- cbind(rnorm(n), rnorm(n)) # 50, 2 
    
    # B from equation 19 
    bsX.true <- cbind(rnorm(kX), rnorm(kX)) # 100, 2
    bX.true <- cbind(rnorm(kX), rnorm(kX), rnorm(kX)) # 100, 3
    bsy.true <- cbind(rnorm(ky), rnorm(ky)) # 40, 2
    by.true <- cbind(rnorm(ky), rnorm(ky)) # 40, 2
    
    dw.z <- diag(c(2, 1)) # 2 by 2 diagonal matrix with 2, 1 on diag 
    dw.y <- diag(c(4, 2)) # 2 by 2 diagonal matrix with 4, 2 on diag --> D(1) = (4,2)
    dw.X <- diag(c(4, 2, 2)) # 3 by 3 diagonal matrix with 4,2,2 on diag --> D(2) = (4,2,2)
    
    # ------------------------------------------------------------------------------------------------------------
    # Equation 19 Y(2)
    # 50 by 100 matrix 
    # ------------------------------------------------------------------------------------------------------------
    X1 <- zs.true %*% dw.z %*% t(bsX.true) +  zX.true %*% dw.X %*% t(bX.true) + 2 * matrix(rnorm(n * kX), nr = n) 
    
    # ------------------------------------------------------------------------------------------------------------
    # Equation 18 Y(1)
    # 50 by 40 matrix 
    # ------------------------------------------------------------------------------------------------------------
    y1 <- zs.true %*% dw.z %*% t(bsy.true) +  zy.true %*% dw.y %*% t(by.true) + 2 * matrix(rnorm(n * ky), nr = n)
    

    # ------------------------------------------------------------------------------------------------------------
    # Notice here b1 is a MD2S obejct
    # ------------------------------------------------------------------------------------------------------------
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
    
    nt <- detectCores()  # 8 core for my poor little laptop
    # from parallel package 
    # Attempt to detect the number of CPU cores on the current host.
    cl <- makeCluster(nt)
    # socket cluster with 8 nodes on host ‘localhost’
    registerDoParallel(cl)
    # register the parallel backend with the foreach package.
    
    # nperm = 200 
    results <- foreach(i = 1:nperm, .export = c('MD2S', 'ginv', 'irlba')) %dopar% {  
      
      # ------------------------------------------------------------------------------------------------------------
      # X1.perm is generated by sampling X1 by column without replacement of the same length
      # So like col 1 in X1.perm should have the exact same numbers as col 1 in X1, but in random order 
      # if we run sort(X1[,1]) == sort(X1.perm[,1]),
      # it should only return a bunch of TRUEs 
      # 
      # sample.mat() is user-defined in MD2S.R
      # Margin = 2 is by columns 
      # apply(X,2,FUN=function(x) sample(x,length(x),replace=FALSE))
      # ------------------------------------------------------------------------------------------------------------
      X1.perm <- sample.mat(sample.mat(X1)) # 50 by 100 matrix

      # y1.perm is generated similarly based on y1 
      y1.perm <- sample.mat(sample.mat(y1)) # 50 by 40 matrix
      
      # MD2S object based on the permuted data 
      b1.perm <- MD2S(X = X1.perm, y = y1.perm, dim = 5)
      
      # ------------------------------------------------------------------------------------------------------------
      # From MD2S object above, we extract the lz, lz.X, lz.y 
      # ------------------------------------------------------------------------------------------------------------
      # output list of transposed matrices 
      out <- list("lz" = t(as.matrix(b1.perm$lz)), 
                  "lz.X" = t(as.matrix(b1.perm$lz.X)),
                  "lz.y" = t(as.matrix(b1.perm$lz.y))
      )
      out
    }
    
    
    stopCluster(cl)
    
    # ------------------------------------------------------------------------------------------------------------
    # initiate lists to store outputs 
    # ------------------------------------------------------------------------------------------------------------
    lz.run <- lz.X.run <- lz.y.run <- list()
    
    # ------------------------------------------------------------------------------------------------------------
    # length of results == nperm (= 200 here)
    # "results" consists of a list of lists containing MD2S results ($lz, $lz.X, $lz.y)
    # All items are of length 5, like results[[1]]$lz has 5 elements for the 5 dimensions 
    # ------------------------------------------------------------------------------------------------------------
    for(i in 1:length(results)) {
      # lz.run contains all lz values from results 
      lz.run[[i]] <- results[[i]]$lz
      # lz.X.run contains all lz.X values from results 
      lz.X.run[[i]] <- results[[i]]$lz.X
      # lz.y.run contains all lz.y values from results 
      lz.y.run[[i]] <- results[[i]]$lz.y
    }
    
    
    # ------------------------------------------------------------------------------------------------------------
    # Cleaning up: 
    # 3 lines below convert the list of lists to matrices of dimension 200 X 5
    # ------------------------------------------------------------------------------------------------------------
    lz.run2  <- do.call('rbind', lz.run)
    lz.X.run2 <- do.call('rbind', lz.X.run)
    lz.y.run2 <- do.call('rbind', lz.y.run)
    
    
    # pz1<-pz2<-pz2.y<-pz2.X<-NULL
    pz2<-pz2.y<-pz2.X<-NULL
    dims <- 5
    
    # ------------------------------------------------------------------------------------------------------------
    # check p values 
    # b1 is the TRUE simulated dataset --> used here to check p value for .1 threshold 
    # b1.perm is the permuted ones 
    # ------------------------------------------------------------------------------------------------------------
    for(i in 1:dims){
      # mean of a list of booleans from the inequality 
      pz2[i]<- mean(b1$lz[i] < lz.run2[,i]) 
      pz2.y[i] <- mean(b1$lz.y[i] < lz.y.run2[,i])
      pz2.X[i] <- mean(b1$lz.X[i] < lz.X.run2[,i])
    }
    
    # p values 
    pvals.s <- pz2
    pvals.X <- pz2.X
    pvals.y <- pz2.y
    
    # pvals.s 
    # pvals.X
    # pvals.y
    
    scenario2 <- 1
    # bind results to a list 
    results.curr<-c(n,kX,scenario2,pvals.s,pvals.X,pvals.y)
    # update the names for the list 
    names(results.curr)<-c("N","kX","scen2",paste("ps_",1:5,sep=""),paste("pX_",1:5,sep=""),paste("py_",1:5,sep=""))
    
    results.all <- results.curr
    
    # output path 
    # file names are generated randomly using uniform distribution
    # My concern is that when we have lots of iterations, there might be duplicated file names 
    # Could we do incrementing file names using a counter? 
    ran.name <- paste0("./results_Panel_A/output_",round(runif(1),10)*1e10, "_Panel_A.RData")
    # save the results from current iteration
    save(results.all, file = ran.name)
  }
}







## ---------------------------------------------------------------------------------------------------
## Aggregate results:
## The code below basically combines the intermediate .RData files into one aggregate file
## I also changed the dir path to simplify things here 
## ---------------------------------------------------------------------------------------------------

files.all<-list.files('./results_Panel_A/', pattern = "output_")
if(!exists("out.all")){
  out.all<-NULL
  for(i in files.all) {
    load(paste('./results_Panel_A/', i, sep = ""))
    out.all<-rbind(out.all,results.all)
  }
}

# save files as one
save(out.all, file = './results_Panel_A/outall_Panel_A.RData')

# # remove intermediate files
# for(i in files.all) {
#   file.remove(paste("./results_Panel_A/", i, sep = ""))
# }
## ---------------------------------------------------------------------------------------------------
