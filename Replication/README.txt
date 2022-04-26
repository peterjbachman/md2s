## --------------------------------------------------------------------------------
## 	Replication Files for: 
##	"Scaling Data from Multiple Sources"
##	Authors: Ted Enamorado, Gabriel Lopez-Moctezuma, Marc Ratkovic
##
##      Abstract:
##      We introduce a method for scaling two data sets from different sources. The
##      proposed method estimates a latent factor common to both datasets as
##      well as an idiosyncratic factor unique to each. In addition, it offers a
##      flexible modeling strategy that permits the scaled locations
##      to be a function of covariates, and efficient
##      implementation allows for inference through resampling.
##      A simulation study shows that our proposed method
##      improves over existing alternatives in capturing the variation common
##      to both datasets, as well as the latent factors specific to each.
##      We apply our proposed method to vote and speech data from
##      the 112th U.S. Senate. We recover a shared subspace that
##      aligns with a standard ideological dimension running from
##      liberals to conservatives, while recovering the words most associated with
##      each senator's location.  In addition, we estimate a word-specific
##      subspace that ranges from national security to budget concerns, and a
##      vote-specific subspace with Tea Party senators on one extreme
##      and senior committee leaders on the other.
##      
## --------------------------------------------------------------------------------

## --------------------------------------------------------------------------------
## 	General Information:
## --------------------------------------------------------------------------------
##
## 	1. Computational Requirements: 
##
##	   There are no specific computational requirements to run the code. 
##	   Any standard personal computer (multi-core) with R installed will 
##         be able to reproduce the results presented in the paper. However,
##	   since we conducted a large number of simulation exercises, we used 
##         Della -- a computer cluster at Princeton University (for more information see:
##         https://researchcomputing.princeton.edu/systems-and-services/available-systems/della).
##         We used 20 cores per node (with 128GB of RAM) and it took around 
##         50 hours to run the code as indicated in the paper.
##
## 	   At Code Ocean, we used AWS EC2 instance with 16 cores with a total 
##         50GB or disk space and 120GB of RAM. R version: 3.6.0.
##
## 	2. Running Time Code Ocean: 
##	   
##	3. List of Files containing the Figures:
##
##		Main Text:
##
##			  - Figure1.pdf         (Figure 1)
##			  - Figure2_Left.pdf    (Figure 2 left panel)
##			  - Figure2_Right.pdf   (Figure 2 right panel)
##			  - Figure3_Panel_A.pdf (Figure 3 panel (a))
##			  - Figure3_Panel_B.pdf (Figure 3 panel (b))
##			  - Figure3_Panel_C.pdf (Figure 3 panel (c))
##			  - Figure4.pdf         (Figure 4)
##			  - Figure5.pdf         (Figure 5)
##			  - Figure6.pdf         (Figure 6)
##
##		Supplemental Appendix: 
##
##			  - SA_Figure_1A.pdf (Figure 1 panel (a))##			  - SA_Figure_1B.pdf (Figure 1 panel (b))##			  - SA_Figure_2A.pdf (Figure 2 panel (a))##			  - SA_Figure_2B.pdf (Figure 1 panel (b))##			  - SA_Figure_2C.pdf (Figure 1 panel (c))##			  - SA_Figure_3A.pdf (Figure 3 panel (a))##			  - SA_Figure_3B.pdf (Figure 3 panel (b))## 			  - SA_Figure_3C.pdf (Figure 3 panel (c))##			  - SA_Figure_4.pdf  (Figure 4)##			  - SA_Figure_5.pdf  (Figure 5)##			  - SA_Figure_6.pdf  (Figure 6)##			  - SA_Figure_7A.pdf (Figure 7 panel (a))## 			  - SA_Figure_7B.pdf (Figure 7 panel (b))##			  - SA_Figure_7C.pdf (Figure 7 panel (c))##			  - SA_Figure_8.pdf  (Figure 8)##			  - SA_Figure_9.pdf  (Figure 9)##			  - SA_Figure_10.pdf (Figure 10)
##
##	4. List of Files containing the Tables:
##
##		Main Text: 
##
##			  - Table_1.txt (Table 1)####		Supplemental Appendix: 
##
##			  - SA_Table_1.txt (Table 1)##			  - SA_Table_2.txt (Table 2)##			  - SA_Table_3.txt (Table 3)##			  - SA_Table_4.txt (Table 4)##			  - SA_Table_5.txt (Table 5)##			  - SA_Table_6.txt (Table 6)##			  - SA_Table_7.txt (Table 7)##
## 	3. Acknowledgments:
##
##	   We are thankful to Simon Heuberger (from Political Analysis) 
##         and to Shahar Zaks (from Code Ocean) for all their guidance
##         with setting the Code Ocean capsule containing our replication 
##	   materials.
##
##      Note 1: The file Run.R call all the necessary files to reproduce every 
##              result in the main text of the paper and in the Supplemental Appendix
##
##      Note 2: 00_MD2S.R in the root folder contains the necessary functions 
##              to implement our proposed method.
##
##      Note 3: To reproduce the code at Code Ocean, we reduced the number of	
##		simulations, permutations, and bootstrap samples from 1000 to
##		200. 
##              
##      Note 4: The tables produced from the code are saved as .txt documents
##		that can be open with any text editor. To format the tables
##		as they appear in the paper, we have used latex and these .txt
##		files as the main input.
##              
## --------------------------------------------------------------------------------

