# md2s


The file gitcommands.txt contains the git commands needed to make changes or checkout your own branch on github.

The instructions.txt file just contains what Ted sent us to work on.

The Replication folder here contains ALL files that we need to work on for the MD2S function and the permutation test.
We don't need to have everything up on the repo for this project. 


PermutationTest.R
This file contains the fully commented version of 01_Create_Data_Panel_A.R for the permutation test. 
Users should be able to specify parameters like ky, n, kX, dims, to replicate Panel A, B, or C, and use for their own application. 

MD2S_original.R
This is just a copy of 00_MD2S.R file. This is the file we need to focus on when building helpfiles for MD2S. Make sure to understand what the MD2S object looks like, what parameters are needed as input.

02_Figure_3.R
This is a copy of plot_Fig3.R file, which plots Panel A in Figure 3. 
This is used just for checking whether our results makes sense, and whether we can replicate Ted's results.
Eventually, we want to wrap this into a function that allows users to plot from the R package. 

Run.R 
Basically sources 00_MD2S.R and runs 01_Create_Data_Panel_A.R
We don't really need it at the current stage. We can just run the code inside each R file.






