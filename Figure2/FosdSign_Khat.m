% this script creates the Figure 2 (left panel)


clear all
close all
clc

n=10000;
lambdag1=0.10; %baseline network parameter
fosdsize=0.005:0.025:0.3; % size of the fosd
lambdagfosd=lambdag1+fosdsize;
betaa=[0.2];
khat=[10:10:2000]; % khat
khat=[1 khat];


for j=1:length(khat)
    for i=1:length(lambdagfosd)
        clc
        display(['lambda ' num2str(lambdagfosd(i)) ' - khat ' num2str(khat(j))]);
        [maxpi_orig(j,i),locmax_orig(j,i),maxL_orig(j,i),maxGAMMA_orig(j,i),maxphi_orig(j,i),maxb_orig(j,i),maxpi_FOSD(j,i),maxL_FOSD(j,i),maxGAMMA_FOSD(j,i),...
            maxphi_FOSD(j,i),maxb_FOSD(j,i),is_FOSD_Larger(j,i),Welfare_orig(j,i),Welfare_FOSD(j,i),diff_avg_orig(j,i),diff_avg_FOSD(j,i)]=CreateNtwFosd_ComputeMax(n,lambdag1,lambdagfosd(i),betaa,khat(j));   
    end
    save('ExplorationF_Rnd_G_Uniform_sizeFOSD_khat2.mat');
end
