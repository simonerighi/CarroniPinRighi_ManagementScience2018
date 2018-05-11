% Creates the Figure 2 - Right panel, where an exemplificatory  baseline
% and FOSD PDFs is depicted.
% To visualize the output:
% 1. Open CreateNtwFosd_ComputeMax
% 2. Set printpdfs=1 (line 16).


clear all
close all


n=10000;
lambdag1=0.1;
fosdsize=0.1;
lambdagfosd=lambdag1+fosdsize;
betaa=0.2;
khat=250;

CreateNtwFosd_ComputeMax(n,lambdag1,lambdagfosd,betaa,khat);