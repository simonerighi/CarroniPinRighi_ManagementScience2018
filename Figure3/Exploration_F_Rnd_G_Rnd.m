% Creates the Data for the Figure 3, for the figure in which both f(k) and g(k) are random networks.

clear all
close all
clc

n=10000;
lambdag1=[0.05:0.01:0.95];
lambdagfosd=lambdag1+0.025;
betaa=[0.5:0.01:0.95];

matlabpool open
for i=1:length(lambdag1)
    parfor j=1:length(betaa)
        clc
        display(['gamma ' num2str(lambdag1(i))]);
        display(['beta ' num2str(betaa(j))]);
        [maxpi_orig(i,j),locmax_orig(i,j),maxL_orig(i,j),maxGAMMA_orig(i,j),maxphi_orig(i,j),maxb_orig(i,j),maxpi_FOSD(i,j),maxL_FOSD(i,j),maxGAMMA_FOSD(i,j),...
            maxphi_FOSD(i,j),maxb_FOSD(i,j),is_FOSD_Larger(i,j),Welfare_orig(i,j),Welfare_FOSD(i,j),diff_avg_orig(i,j),diff_avg_FOSD(i,j)]=fct_F_Rnd_G_Rnd(n,lambdag1(i),lambdagfosd(i),betaa(j));   
    end
    save('ExplorationF_Rnd_G_Rnd.mat');
end
matlabpool close

clear all
load ExplorationF_Rnd_G_Rnd.mat

diff=maxpi_FOSD-maxpi_orig;

figure(1)
contour(betaa,lambdag1,diff,'ShowText','on');
set(gca,'YDir','normal')
ylabel('\lambda_g','FontSize',18);
xlabel('\beta','FontSize',18);
colorbar
title({'Outdegree:Rand - Indegree:Rand' '\pi^*_{FOSD}-\pi^*_{Orig}'},'Fontsize',20)
print -depsc Explore_FisRand_GisRand.eps


diff_sign=(diff>0).*2-1;

figure(2)
contour(betaa,lambdag1,diff_sign,'ShowText','on');
set(gca,'YDir','normal')
ylabel('\lambda_g','FontSize',18);
xlabel('\beta','FontSize',18);
colorbar
title({'Outdegree:Rand - Indegree:Rand';'Sign of \pi^*_{FOSD}-\pi^*_{Orig}'},'Fontsize',20)
print -depsc Explore_FisRand_GisRand_Sign.eps


diff_Welfare=Welfare_FOSD-Welfare_orig;

figure(3)
contour(betaa,lambdag1,diff_Welfare,'ShowText','on');
set(gca,'YDir','normal')
ylabel('\lambda_g','FontSize',18);
xlabel('\beta','FontSize',18);
colorbar
title({'Outdegree:Rand - Indegree:Rand' 'W^*_{FOSD}-W^*_{Orig}'},'Fontsize',20)
print -depsc Explore_Welfare_FisRand_GisRand.eps


diff_Welfare_sign=(diff_Welfare>0).*2-1;

figure(4)
contour(betaa,lambdag1,diff_Welfare_sign,'ShowText','on');
set(gca,'YDir','normal')
ylabel('\lambda_g','FontSize',18);
xlabel('\beta','FontSize',18);
colorbar
title({'Outdegree:Rand - Indegree:Rand';'Sign of W^*_{FOSD}-W^*_{Orig}'},'Fontsize',20)
print -depsc Explore_Welfare_FisRnd_GisRand_Sign.eps
