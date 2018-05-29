% This generates results of Figure 4
% Carroni Pin Righi (2018), submitted for pubblication in Management
% Sciences


clear all
close all
clc
n=10000;
gammag1=2.0:0.01:3;
betaa=0.05;

maxpi_orig=zeros(length(betaa),length(gammag1));
maxL_orig=zeros(length(betaa),length(gammag1));
maxGAMMA_orig=zeros(length(betaa),length(gammag1));
maxphi_orig=zeros(length(betaa),length(gammag1));
klowerbest=zeros(length(betaa),length(gammag1));
klowerbestrelative=zeros(length(betaa),length(gammag1));


for i=1:length(betaa)
    for j=1:length(gammag1)
       display(['beta ' num2str(betaa(i)) ' - Gamma ' num2str(gammag1(j))])
       [maxpi_orig(i,j),maxL_orig(i,j),maxGAMMA_orig(i,j),maxphi_orig(i,j),klowerbest(i,j),klowerbestrelative(i,j),klowerbest_only_public(i,j)]=Create_Network_ComputeMaximum(n,gammag1(j),betaa(i));
      
    end
end

save('FSF_GRND_mixed_changeF')


% run only the code below to generate figure from available data, if you
% run everything it will re-create the data for the graph (depending on
% input it may take some time).
clear all
close all
clc
load FSF_GRND_mixed_changeF

figure(1)
semilogy(gammag1,klowerbest(1,:),'b-');
hold on
semilogy(gammag1,klowerbest_only_public(1,:),'r-');
legend('Public and Private WOM','Only Public WOM')
title({'Optimal $$\underline{k}$$ - Different $$\gamma_f$$';'Indegree:Random - Outdegree:SF'},'Interpreter','latex','FontSize',18)
xlabel('\gamma_f','FontSize',18)
ylabel('$$\underline{k}^*$$','Interpreter','latex','FontSize',18)
ylim([0 100])
print -depsc MixedChangeF_Grand_FSF.eps
