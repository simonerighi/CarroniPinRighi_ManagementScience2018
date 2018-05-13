% This script computes the effect of a MPS on f(k) AND a MPS on g(k), when both g(k) and f(k)
% are random networks. The beta is adapted to fulfill the consistency
% condition between the out and the in- degree distribution.
% study of the mean preserving spread of f
% the effect of a fosd for different lambda_f and lambda_g is computed and
% made into a graph
clear all
close all
clc



n=100000;
plink_f=0.05:0.05:0.95; % lambda_f
plink_g=0.05:0.05:0.95; % lambda_g
p1=1/2;
p2=1/2;
cost=0.06;
kmin=1;
kmax=n-2;
delta_gradi=30; % size of the MPS in degrees

warning off

effect_f=zeros(length(plink_f),length(plink_g));
klower_star_mps_f=zeros(length(plink_f),length(plink_g));
L_star_mps_f=zeros(length(plink_f),length(plink_g));
GAMMA_star_mps_f=zeros(length(plink_f),length(plink_g));
phi_star_mps_f=zeros(length(plink_f),length(plink_g));
b_star_mps_f=zeros(length(plink_f),length(plink_g));
klower_star_orig_f=zeros(length(plink_f),length(plink_g));
L_star_orig_f=zeros(length(plink_f),length(plink_g));
GAMMA_star_orig_f=zeros(length(plink_f),length(plink_g));
phi_star_orig_f=zeros(length(plink_f),length(plink_g));
b_star_orig_f=zeros(length(plink_f),length(plink_g));
profit_mps_f=zeros(length(plink_f),length(plink_g));
profit_orig_f=zeros(length(plink_f),length(plink_g));

effect_g=zeros(length(plink_f),length(plink_g));
klower_star_mps_g=zeros(length(plink_f),length(plink_g));
L_star_mps_g=zeros(length(plink_f),length(plink_g));
GAMMA_star_mps_g=zeros(length(plink_f),length(plink_g));
phi_star_mps_g=zeros(length(plink_f),length(plink_g));
b_star_mps_g=zeros(length(plink_f),length(plink_g));
klower_star_orig_g=zeros(length(plink_f),length(plink_g));
L_star_orig_g=zeros(length(plink_f),length(plink_g));
GAMMA_star_orig_g=zeros(length(plink_f),length(plink_g));
phi_star_orig_g=zeros(length(plink_f),length(plink_g));
b_star_orig_g=zeros(length(plink_f),length(plink_g));
profit_mps_g=zeros(length(plink_f),length(plink_g));
profit_orig_g=zeros(length(plink_f),length(plink_g));

diff_observed_f=zeros(length(plink_f),length(plink_g))-inf;
diff_observed_g=zeros(length(plink_f),length(plink_g))-inf;

betaa=zeros(length(plink_f),length(plink_g));

for i=1:length(plink_f) % for all lambda_f
    for j=1:length(plink_g) % for all lambda_g
         gf=plink_f(i);
         gg=plink_g(j);
         display(['\gamma f= ' num2str(gf)]);
         display(['\gamma g= ' num2str(gg)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate MPS network of f
        kval=[1:1:n-2];
        fdik=binopdf(kval,n-2,gf);
        meank_f=sum(fdik.*kval);
        mean_f_mps_mtx(i,j)=meank_f;
        
        % remove unnecessary 0s in the extremes of the distribution.
        dummyfdik=fdik;
        if find(dummyfdik==0,1)==1; kmin_mps_f=find(dummyfdik>0,1); dummyfdik(1:kmin_mps_f-1)=dummyfdik(1:kmin_mps_f-1)+1; else; kmin_mps_f=1; end;
        if dummyfdik(end)==0; kmax_mps_f=find(dummyfdik==0,1)-1; else; kmax_mps_f=n-2; end
        clear dummyfdik
        fdik_mps=fdik(kmin_mps_f:kmax_mps_f);
        kval_mps_f=kmin_mps_f:1:kmax_mps_f;

        clear fdik
        clear kval

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generated MPS network for g
        kval=[1:1:n-2];
        gdik=binopdf(kval,n-2,gg);
        meank_g=sum(gdik.*kval);
        mean_g_mps_mtx(i,j)=meank_g;

         % remove unnecessary 0s in the extremes of the distribution.
        dummygdik=gdik;
        if find(dummygdik==0,1)==1; kmin_mps_g=find(dummygdik>0,1); dummygdik(1:kmin_mps_g-1)=dummygdik(1:kmin_mps_g-1)+1; else; kmin_mps_g=1; end;
        if dummygdik(end)==0; kmax_mps_g=find(dummygdik==0,1)-1; else; kmax_mps_g=n-2; end; %length(dummyfdik)
        clear dummyfdik
        gdik_mps=gdik(kmin_mps_g:kmax_mps_g); % mi salvo k e fdik per il mps
        kval_mps_g=kmin_mps_g:1:kmax_mps_g;
        clear simtillN 
        clear gdik
        clear kval

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generated baseline network f(k) which is more concentrated of delta_degrees
        avgdegree=meank_f;
        %Take f(k) values from the original network and collect all the mass that corresponds to k,s.t. (lambda*n-delta) > k > (lambda*n+delta)
        fdik_orig_f=fdik_mps;
        kmin_orig_f=round(avgdegree-delta_gradi);
        kmax_orig_f=round(avgdegree+delta_gradi);
        kval_orig_f=kmin_orig_f:1:kmax_orig_f;
        Garthered_mass=sum(fdik_orig_f(kval_mps_f<kmin_orig_f))+sum(fdik_orig_f(kval_mps_f>kmax_orig_f));
        fdik_orig_f(kval_mps_f<kmin_orig_f)=0;
        fdik_orig_f(kval_mps_f>kmax_orig_f)=0;

        to_redistribute=Garthered_mass/length(kval_orig_f); % Redistribute this same mass equally divided among each k value with k,s.t. (lambda*n-delta) < k < (lambda*n+delta).
        fdik_orig2_f=fdik_orig_f(fdik_orig_f>0);


        % i redistribute the collected mass proportionally to the mass
        % remaining.
        fdik_orig3_f=fdik_orig2_f+to_redistribute;
        clear fdik_orig
        fdik_orig=fdik_orig3_f;
        clear fdik_orig3
        clear fdik_orig2
        mean_f_orig_mtx(i,j)=sum(fdik_orig.*kval_orig_f);
        clear to_redistribute
        clear Garthered_mass

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate baseline network g(k), also more concentrated of
        % delta_gradi degrees.
        
        avgdegree=meank_g;
        %Take g(k) values from the original network and collect all the mass that corresponds to k,s.t. (lambda*n-delta) > k > (lambda*n+delta)
        gdik_orig_g=gdik_mps;
        kmin_orig_g=round(avgdegree-delta_gradi);
        kmax_orig_g=round(avgdegree+delta_gradi);
        kval_orig_g=kmin_orig_g:1:kmax_orig_g;
        Garthered_mass=sum(gdik_orig_g(kval_mps_g<kmin_orig_g))+sum(gdik_orig_g(kval_mps_g>kmax_orig_g));
        gdik_orig_g(kval_mps_g<kmin_orig_g)=0;
        gdik_orig_g(kval_mps_g>kmax_orig_g)=0;

        to_redistribute=Garthered_mass/length(kval_orig_g); % Redistribute this same mass equally divided among each k value with k,s.t. (lambda*n-delta) < k < (lambda*n+delta).
        gdik_orig2_g=gdik_orig_g(gdik_orig_g>0);


        % i redistribute the collected mass proportionally to the mass
        % remaining.
        gdik_orig3_g=gdik_orig2_g+to_redistribute;
        clear gdik_orig_g
        gdik_orig=gdik_orig3_g;
        clear gdik_orig3_g
        clear gdik_orig2_g
        
        mean_g_orig_mtx(i,j)=sum(gdik_orig.*kval_orig_g);
       
        
        % compute the betaa such that the consistency condition is
        % fulfilled.
        betaa(i,j)=meank_f/(meank_g+meank_f);
        
        diff_observed_f=mean_f_mps_mtx(i,j)-mean_f_orig_mtx(i,j);
        diff_observed_g=mean_g_mps_mtx(i,j)-mean_g_orig_mtx(i,j);

        % compute effect on profit of mps on F
        [effect_f(i,j),klower_star_mps_f(i,j),L_star_mps_f(i,j),GAMMA_star_mps_f(i,j),phi_star_mps_f(i,j),b_star_mps_f(i,j),klower_star_orig_f(i,j),L_star_orig_f(i,j),GAMMA_star_orig_f(i,j),...
            phi_star_orig_f(i,j),b_star_orig_f(i,j),profit_mps_f(i,j),profit_orig_f(i,j)]=...
            effect_of_mps_f(kmax_mps_f,kmin_mps_f,kmax_orig_f,kmin_orig_f,kval_mps_f,kval_orig_f,kval_orig_g,gdik_orig,fdik_orig,fdik_mps,...
            cost,betaa(i,j),p1,p2);

        % compute effect on profit of mps on G
        [effect_g(i,j),klower_star_mps_g(i,j),L_star_mps_g(i,j),GAMMA_star_mps_g(i,j),phi_star_mps_g(i,j),b_star_mps_g(i,j),klower_star_orig_g(i,j),L_star_orig_g(i,j),GAMMA_star_orig_g(i,j),...
            phi_star_orig_g(i,j),b_star_orig_g(i,j),profit_mps_g(i,j),profit_orig_g(i,j)]=effect_of_mps_g(...
            kmax_mps_g,kmin_mps_g,kmax_orig_f,kmin_orig_f,kval_mps_g,kval_orig_g,kval_orig_f,gdik_orig,fdik_orig,gdik_mps,...
            cost,betaa(i,j),p1,p2,fdik_mps); 

    end
end

% print differences in profits for MPS in g
figure(13)
imagesc(plink_f,plink_g,profit_mps_g-profit_orig_g)
set(gca,'YDir','normal')
xlabel('\lambda_f','FontSize',18);
ylabel('\lambda_g','FontSize',18);
colorbar
title('Effect of a MPS on g(k) on Difference in \pi*','Fontsize',20)
print -depsc Random_FdifferentG_pi_G.eps


% print differences in profits for MPS in f
figure(14)
imagesc(plink_f,plink_g,profit_mps_f-profit_orig_f)
set(gca,'YDir','normal')
xlabel('\lambda_f','FontSize',18);
ylabel('\lambda_g','FontSize',18);
colorbar
title('Effect of a MPS on f(k) on Difference in \pi*','Fontsize',20)
print -depsc Random_FdifferentG_pi_F.eps

    
% figure(1)
% imagesc(plink_f,plink_g,effect_f)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Sign of the effect of a MPS on f(k) on profits','Fontsize',20)
% print -depsc Random_FdifferentG_effect_F.eps
% 
% 
% figure(2)
% imagesc(plink_f,plink_g,klower_star_mps_f-klower_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on K*','Fontsize',20)
% print -depsc Random_FdifferentG_K_F.eps
% 
% 
% figure(3)
% imagesc(plink_f,plink_g,L_star_mps_f-L_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on L*','Fontsize',20)
% print -depsc Random_FdifferentG_L_F.eps
% 
% 
% figure(4)
% imagesc(plink_f,plink_g,GAMMA_star_mps_f-GAMMA_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on \Gamma*','Fontsize',20)
% print -depsc Random_FdifferentG_GAMMA_F.eps
% 
% 
% figure(5)
% imagesc(plink_f,plink_g,phi_star_mps_f-phi_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on \phi*','Fontsize',20)
% print -depsc Random_FdifferentG_phi_F.eps
% 
% figure(6)
% imagesc(plink_f,plink_g,b_star_mps_f-b_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on b*','Fontsize',20)
% print -depsc Random_FdifferentG_b_F.eps
% 
% 
% figure(7)
% imagesc(plink_f,plink_g,effect_g)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Sign of the effect of a MPS on g(k) on profits','Fontsize',20)
% print -depsc Random_FdifferentG_effect_G.eps
% 
% 
% figure(8)
% imagesc(plink_f,plink_g,klower_star_mps_g-klower_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on K*','Fontsize',20)
% print -depsc Random_FdifferentG_K_G.eps
% 
% figure(9)
% imagesc(plink_f,plink_g,L_star_mps_g-L_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on L*','Fontsize',20)
% print -depsc Random_FdifferentG_L_G.eps
% 
% 
% figure(10)
% imagesc(plink_f,plink_g,GAMMA_star_mps_g-GAMMA_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on \Gamma*','Fontsize',20)
% print -depsc Random_FdifferentG_GAMMA_G.eps
% 
% 
% figure(11)
% imagesc(plink_f,plink_g,phi_star_mps_g-phi_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on \phi*','Fontsize',20)
% print -depsc Random_FdifferentG_phi_G.eps
% 
% 
% figure(12)
% imagesc(plink_f,plink_g,b_star_mps_g-b_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on b*','Fontsize',20)
% print -depsc Random_FdifferentG_b_G.eps
%
% 
% figure(15)
% imagesc(plink_f,plink_g,diff_observed_f)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Numerical Difference mean(ORIG_f)-mean(MPS_f)','Fontsize',20)
% print -depsc Random_FequalG_DiffObs_f.eps
% 
% figure(16)
% imagesc(plink_f,plink_g,diff_observed_g)
% set(gca,'YDir','normal')
% xlabel('\lambda_f','FontSize',18);
% ylabel('\lambda_g','FontSize',18);
% colorbar
% title('Numerical Difference mean(ORIG_g)-mean(MPS_g)','Fontsize',20)
% print -depsc Random_FequalG_DiffObs_g.eps
