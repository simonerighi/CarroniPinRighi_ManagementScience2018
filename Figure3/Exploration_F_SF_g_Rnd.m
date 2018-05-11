% Creates the Data for the Figure 3, for the figure in which f(k) is a scale-free network while g(k) is a random network

clear all
close all
clc

n=10000;
gammag1=[0.05:0.01:0.5]; % baseline gamma
gammagfosd=gammag1+0.025; 
betaa=[0.5:0.01:0.95]; %beta


for i=1:length(gammag1)
    for j=1:length(betaa)
        clc
        display(['gamma ' num2str(gammag1(i))]);
        display(['beta ' num2str(betaa(j))]);
        [maxpi_orig(i,j),locmax_orig(i,j),maxL_orig(i,j),maxGAMMA_orig(i,j),maxphi_orig(i,j),maxb_orig(i,j),maxpi_FOSD(i,j),maxL_FOSD(i,j),...
            maxGAMMA_FOSD(i,j),...
            maxphi_FOSD(i,j),maxb_FOSD(i,j),is_FOSD_Larger(i,j),Welfare_orig(i,j),Welfare_FOSD(i,j),diff_avg_orig(i,j),diff_avg_FOSD(i,j)]...
            =fct_F_sf_G_Rnd(n,gammag1(i),gammagfosd(i),betaa(j));
    end
    save('ExplorationF_SF_G_Rnd.mat');
end

