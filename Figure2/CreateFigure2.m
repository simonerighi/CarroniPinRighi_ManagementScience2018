% this script creates the Figure 2 (left panel) once the numerical
% simulations are run launching the script FosdSign_Khat.m

clear all
close all
load ExplorationF_Rnd_G_Uniform_sizeFOSD_khat2.mat

diff=maxpi_FOSD-maxpi_orig;

sign_diff=zeros(size(diff));

for i=1:size(diff,1)
    for j=1:size(diff,2)
        if diff(i,j)>0
            sign_diff(i,j)=1;
        elseif diff(i,j)<0
            sign_diff(i,j)=-1;
        elseif diff(i,j)==0
            sign_diff(i,j)=0;
        else
            sign_diff(i,j)=99;
        end   
    end
end

figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
vctx=[10:10:2000];
vctx=vctx./n;
vcty=[0.005 0.03 0.055 0.08 0.105 0.13 0.155 0.18 0.205 0.23 0.255 0.28];
% Create image
image(vcty,vctx,sign_diff,'Parent',axes1,'CDataMapping','scaled');
%imagesc(khat,fosdsize,sign_diff);
set(gca,'YDir','normal')
xlabel('\Delta \lambda_g','FontSize',18);
ylabel('$$\hat{k}/k_{max}$$','Interpreter','Latex','FontSize',18);%,
%Uncomment the following line to preserve the X-limits of the axes
ylim(axes1,[0.001 0.201]);
%Uncomment the following line to preserve the Y-limits of the axes
xlim(axes1,[0.003 0.293]);

colorbar('peer',axes1,'Position',...
    [0.933799342105263 0.126272912167915 0.0131578947368421 0.750509174689508],...
    'Ticks',[-1 0 1],...
    'TickLabels',{'-1','0','1'},...
    'FontSize',20);

title({'Outdegree:Uniform - Indegree:Rand' '\pi^*_{FOSD}-\pi^*_{Orig}'},'Fontsize',20)

print -depsc FixLambda_FisRand_GisUniform_khat_signdiff_pi.eps

