% Generates the graphs of Figure 3 from available data


clear all
close all
clc

load Exploration_F_Rnd_G_SF.mat

diff_pi=maxpi_FOSD-maxpi_orig;
diff_W=Welfare_FOSD-Welfare_orig;
[sizematr_x,sizematr_y]=size(diff_W);
toprint=zeros(sizematr_x,sizematr_y);
for i=1:sizematr_x
    for j=1:sizematr_y
        if isnan(diff_pi(i,j))==1 || isnan(diff_W(i,j))==1
            toprint(i,j)=0;
        else
            if diff_pi(i,j)>=0 && diff_W(i,j)>=0
                toprint(i,j)=1;
            end
            if diff_pi(i,j)>=0 && diff_W(i,j)<0
               toprint(i,j)=2; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)>=0
               toprint(i,j)=3; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)<0
               toprint(i,j)=4; 
            end   
        end
    end
end

figure(1)
imagesc(betaa,gammag1,toprint); % 'ShowText','on'
set(gca,'YDir','normal')
ylabel('\gamma_g','FontSize',18);
xlabel('\beta','FontSize',18);
mymap = [1 1 1
1 0 0
    0 1 0
    0 0 1
    0 0 0 ];
caxis([0, 4])
colormap(mymap)
colorbar('Location','southoutside','Ticks',[0.4,1.20,2,2.8,3.6],'TickLabels',{'NaN','$\Delta(\pi),\Delta(S)>0$','$\Delta(S)<0<\Delta(\pi)$',...
    '$\Delta(\pi)<0<\Delta(S)$','$\Delta(\pi),\Delta(S)<0$'},'TickLabelInterpreter','latex')
title({'Seller''s incentives' 'Outdegree:Rand - Indegree:SF'},'Fontsize',20)
print -depsc IncentiveSigns_Explore_F_Rand_G_SF2.eps
print -dpdf IncentiveSigns_Explore_F_Rand_G_SF2.pdf

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

clear all
close all
clc

load ExplorationF_SF_G_SF.mat

diff_pi=maxpi_FOSD-maxpi_orig;
diff_W=Welfare_FOSD-Welfare_orig;
[sizematr_x,sizematr_y]=size(diff_W);
toprint=zeros(sizematr_x,sizematr_y);
for i=1:sizematr_x
    for j=1:sizematr_y
        if isnan(diff_pi(i,j))==1 || isnan(diff_W(i,j))==1
            toprint(i,j)=0;
        else
            if diff_pi(i,j)>=0 && diff_W(i,j)>=0
                toprint(i,j)=1;
            end
            if diff_pi(i,j)>=0 && diff_W(i,j)<0
               toprint(i,j)=2; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)>=0
               toprint(i,j)=3; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)<0
               toprint(i,j)=4; 
            end   
        end
    end
end

figure(1)
imagesc(betaa,gammag1,toprint); % 'ShowText','on'
set(gca,'YDir','normal')
ylabel('\gamma_g','FontSize',18);
xlabel('\beta','FontSize',18);
mymap = [1 1 1
1 0 0
    0 1 0
    0 0 1
    0 0 0 ];
caxis([0, 4])
colormap(mymap)
colorbar('Location','southoutside','Ticks',[0.4,1.20,2,2.8,3.6],'TickLabels',{'NaN','$\Delta(\pi),\Delta(S)>0$','$\Delta(S)<0<\Delta(\pi)$',...
    '$\Delta(\pi)<0<\Delta(S)$','$\Delta(\pi),\Delta(S)<0$'},'TickLabelInterpreter','latex')
title({'Seller''s incentives' 'Outdegree:SF - Indegree:SF'},'Fontsize',20)
print -depsc IncentiveSigns_Explore_F_SF_G_SF2.eps
print -dpdf IncentiveSigns_Explore_F_SF_G_SF2.pdf


%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

clear all
close all
clc

load ExplorationF_Rnd_G_Rnd.mat

diff_pi=maxpi_FOSD-maxpi_orig;
diff_W=Welfare_FOSD-Welfare_orig;
[sizematr_x,sizematr_y]=size(diff_W);
toprint=zeros(sizematr_x,sizematr_y);
for i=1:sizematr_x
    for j=1:sizematr_y
        if isnan(diff_pi(i,j))==1 || isnan(diff_W(i,j))==1
            toprint(i,j)=0;
        else
            if diff_pi(i,j)>=0 && diff_W(i,j)>=0
                toprint(i,j)=1;
            end
            if diff_pi(i,j)>=0 && diff_W(i,j)<0
               toprint(i,j)=2; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)>=0
               toprint(i,j)=3; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)<0
               toprint(i,j)=4; 
            end   
        end
    end
end

figure(1)
imagesc(betaa,lambdag1,toprint); % 'ShowText','on'
set(gca,'YDir','normal')
ylabel('\lambda_g','FontSize',18);
xlabel('\beta','FontSize',18);
mymap = [1 1 1
1 0 0
    0 1 0
    0 0 1
    0 0 0 ];
caxis([0, 4])
colormap(mymap)
colorbar('Location','southoutside','Ticks',[0.4,1.20,2,2.8,3.6],'TickLabels',{'NaN','$\Delta(\pi),\Delta(S)>0$','$\Delta(S)<0<\Delta(\pi)$',...
    '$\Delta(\pi)<0<\Delta(S)$','$\Delta(\pi),\Delta(S)<0$'},'TickLabelInterpreter','latex')
title({'Seller''s incentives' 'Outdegree:Rand - Indegree:Rand'},'Fontsize',20)
print -depsc IncentiveSigns_Explore_F_Rnd_G_Rnd2.eps
print -dpdf IncentiveSigns_Explore_F_Rnd_G_Rnd2.pdf


%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

clear all
close all
clc

load ExplorationF_SF_G_Rnd.mat

diff_pi=maxpi_FOSD-maxpi_orig;
diff_W=Welfare_FOSD-Welfare_orig;
[sizematr_x,sizematr_y]=size(diff_W);
toprint=zeros(sizematr_x,sizematr_y);
for i=1:sizematr_x
    for j=1:sizematr_y
        if isnan(diff_pi(i,j))==1 || isnan(diff_W(i,j))==1
            toprint(i,j)=0;
        else
            if diff_pi(i,j)>=0 && diff_W(i,j)>=0
                toprint(i,j)=1;
            end
            if diff_pi(i,j)>=0 && diff_W(i,j)<0
               toprint(i,j)=2; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)>=0
               toprint(i,j)=3; 
            end
            if diff_pi(i,j)<0 && diff_W(i,j)<0
               toprint(i,j)=4; 
            end   
        end
    end
end

figure(1)
imagesc(betaa,gammag1,toprint); % 'ShowText','on'
set(gca,'YDir','normal')
ylabel('\lambda_g','FontSize',18);
xlabel('\beta','FontSize',18);
mymap = [1 1 1
1 0 0
    0 1 0
    0 0 1
    0 0 0 ];
caxis([0, 4])
colormap(mymap)
colorbar('Location','southoutside','Ticks',[0.4,1.20,2,2.8,3.6],'TickLabels',{'NaN','$\Delta(\pi),\Delta(S)>0$','$\Delta(S)<0<\Delta(\pi)$',...
    '$\Delta(\pi)<0<\Delta(S)$','$\Delta(\pi),\Delta(S)<0$'},'TickLabelInterpreter','latex')
title({'Seller''s incentives' 'Outdegree:SF - Indegree:Rand'},'Fontsize',20)
print -depsc IncentiveSigns_Explore_F_SF_G_Rnd2.eps
print -dpdf IncentiveSigns_Explore_F_SF_G_Rnd2.pdf








