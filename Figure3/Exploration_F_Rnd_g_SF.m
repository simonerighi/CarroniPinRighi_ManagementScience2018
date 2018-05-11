% Creates the Data for the Figure 3, for the figure in which f(k) is a random network while g(k) is a Scale-Free network

clear all
close all
clc

n=10000;
gammag1=[2.1:0.01:3];
gammagfosd=gammag1-0.025;
betaa=[0.5:0.01:0.95];

for i=1:length(gammag1)
    for j=1:length(betaa)
        clc
        display(['gamma ' num2str(gammag1(i))]);
        display(['beta ' num2str(betaa(j))]);
        [maxpi_orig(i,j),locmax_orig(i,j),maxL_orig(i,j),maxGAMMA_orig(i,j),maxphi_orig(i,j),maxb_orig(i,j),maxpi_FOSD(i,j),maxL_FOSD(i,j),...
            maxGAMMA_FOSD(i,j),...
            maxphi_FOSD(i,j),maxb_FOSD(i,j),is_FOSD_Larger(i,j),Welfare_orig(i,j),Welfare_FOSD(i,j),diff_avg_orig(i,j),diff_avg_FOSD(i,j)]=fct_F_Rnd_G_SF(n,gammag1(i),gammagfosd(i),...
            betaa(j));
    end
    save('Exploration_F_Rnd_G_SF.mat');
end
