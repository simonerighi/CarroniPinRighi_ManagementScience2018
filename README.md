# CarroniPinRighi_ManagementScience2018
This repository contains the simulations of the article from Elias Carroni, Paolo Pin and Simone Righi, submitted for publication to Management Sciences

- Files contained in each folder allow to reproduce the panels of the figure with the same number. 


Figure 2:
- To create the data for the left panel run "FosdSign_Khat.m". The figure can be generated using "CreateFigure2.m" (note that data are already contained in this folder so the figure can be created directly).
- To create the data and the figure of the right panel, use "Khat_showPDF.m" following the instructions contained at the beginning of the file.

Figure 3:
- To create the data for the left panel run "Figure3.m" or the script called by the latter script. 
- Figure can be generated using "Figure3_CreateFigure.m" (note that data are already contained in this folder so the figures can be created directly).


Figure 4: 
- To create the data and print the figure 4 use "MixedStrategy_OptimalK.m".

MPS: This folder contains script and results to create a Mean Preserving Spread (MPS) on a Network and to assess its results on the profits generated by the firm at the optimum.

The folder contains two subfolders: 
  -the first is denominated "Random_Networks" and contains the script and results concerning Mean Preserving Spreads on Bernoulli Random Networks (Erdos-Renyi)
  - the second is denominated "ScaleFree_Networks" and contains the script and results concerning Mean Preserving Spreads on Scale_Free networks Networks (Erdos-Renyi)
  
  In both cases the result of the MPS of a fixed size is computed for different combinations of network parameter values.
  
  In both case, the figures can be reproduced launching the script "MPS_on_F_and_G.m"
