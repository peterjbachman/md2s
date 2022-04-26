## ---------------------------------------------------------------------------------------------------
##
## 		Master file:
##		"Scaling Data from Multiple Sources"
##		Authors: Ted Enamorado, Gabriel Lopez-Moctezuma, Marc Ratkovic
## 
## ---------------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------------
## Load MD2S Function:
## ---------------------------------------------------------------------------------------------------
source("./00_MD2S.R")

##  To create the simulated data for each plot set `create.sim.data <- TRUE'
##  Warning: if you are using a standard laptop computer creating the necessary 
##  datasets necessary to reproduce Figures 1, 2, and 3, and Table 1
##	might take a lot of time. See readme.txt for running times on a
##	standard laptop. A computer cluster would be an option. 

## ---------------------------------------------------------------------------------------------------
## Code we need to run 
## Number of simulations; permutations; and bootstap samples
library(parallel)
nsims <- nperm <- nboot <- 200
# Be aware! This takes a longgggg time to run on your laptop 
source("./01_Create_Data_Panel_A.R")
## ---------------------------------------------------------------------------------------------------

