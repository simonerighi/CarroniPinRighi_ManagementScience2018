% This script computes the effect of a MPS on f(k) AND a MPS on g(k), when both g(k) and f(k)
% are ScaleFree networks. The beta is adapted to fulfill the consistency
% condition between the out and the in- degree distribution.
% study of the mean preserving spread of f
% the effect of a fosd for different lambda_f and lambda_g is computed and
% made into a graph
clear all
close all
clc



n=100000;
gammf=2.00:0.05:3; % gamma f
gammg=2.00:0.05:3; % gamma g
p1=1/2;
p2=1/2;
cost=0.06;
kmin=10;
kmax=n-1;
howclose=0;


effect_f=zeros(length(gammf),length(gammg));
klower_star_mps_f=zeros(length(gammf),length(gammg));
L_star_mps_f=zeros(length(gammf),length(gammg));
GAMMA_star_mps_f=zeros(length(gammf),length(gammg));
phi_star_mps_f=zeros(length(gammf),length(gammg));
b_star_mps_f=zeros(length(gammf),length(gammg));
klower_star_orig_f=zeros(length(gammf),length(gammg));
L_star_orig_f=zeros(length(gammf),length(gammg));
GAMMA_star_orig_f=zeros(length(gammf),length(gammg));
phi_star_orig_f=zeros(length(gammf),length(gammg));
b_star_orig_f=zeros(length(gammf),length(gammg));
profit_mps_f=zeros(length(gammf),length(gammg));
profit_orig_f=zeros(length(gammf),length(gammg));

effect_g=zeros(length(gammf),length(gammg));
klower_star_mps_g=zeros(length(gammf),length(gammg));
L_star_mps_g=zeros(length(gammf),length(gammg));
GAMMA_star_mps_g=zeros(length(gammf),length(gammg));
phi_star_mps_g=zeros(length(gammf),length(gammg));
b_star_mps_g=zeros(length(gammf),length(gammg));
klower_star_orig_g=zeros(length(gammf),length(gammg));
L_star_orig_g=zeros(length(gammf),length(gammg));
GAMMA_star_orig_g=zeros(length(gammf),length(gammg));
phi_star_orig_g=zeros(length(gammf),length(gammg));
b_star_orig_g=zeros(length(gammf),length(gammg));
profit_mps_g=zeros(length(gammf),length(gammg));
profit_orig_g=zeros(length(gammf),length(gammg));

diff_observed_f=zeros(length(gammf),length(gammg))-inf;
diff_observed_g=zeros(length(gammf),length(gammg))-inf;

flag_f=zeros(length(gammf),length(gammg));
flag_g=zeros(length(gammf),length(gammg));


betaa=zeros(length(gammf),length(gammg));

for i=1:length(gammf) % for all gamma_f
    for j=1:length(gammg) % for all  gamma_g
         gf=gammf(i);
         gg=gammg(j);
         display(['\gamma f= ' num2str(gf)]);
         display(['\gamma g= ' num2str(gg)]);
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate MPS network of f
        count=1;
        fdik=zeros(kmax-kmin,1);
        kval=zeros(kmax-kmin,1);
        kk=kmin:kmax;
        simtillN=sum((1./(kk.^gf)));
        for k=kmin:1:kmax
            fdik(count)=(1/k^(gf))/simtillN;
            kval(count)=k;
            count=count+1;
        end
        meank_f=sum(fdik.*kval);
        fdik_mps=fdik; % mi salvo k e fdik per il mps
        kval_mps_f=kval;
        kmin_mps_f=kmin;
        kmax_mps_f=kmax;


        clear simtillN 
        clear fdik
        clear kval

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate MPS network for g
        count=1;
        gdik=zeros(kmax-kmin,1);
        kval=zeros(kmax-kmin,1);
        kk=kmin:kmax;
        simtillN=sum((1./(kk.^gg)));
        for k=kmin:1:kmax
            gdik(count)=(1/k^(gg))/simtillN;
            kval(count)=k;
            count=count+1;
        end
        meank_gg=sum(gdik.*kval);
        gdik_mps=gdik; % mi salvo k e fdik per il mps
        kval_mps_g=kval;
        kmin_mps_g=kmin;
        kmax_mps_g=kmax;
        
        clear simtillN 
        clear fdik
        clear kval

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generated baseline network f(k) which is more concentrated of delta_degrees
        kmin_round=kmin+1; % to have a more concentrated network i reduce of 1 the span of possible ks
        count=1;
        fdik=zeros(kmax-kmin_round,1);
        kval=zeros(kmax-kmin_round,1);
        kk=kmin_round:kmax;
        simtillN=sum((1./(kk.^gf)));
        for k=kmin:1:kmax
            fdik(count)=(1/k^(gf))/simtillN;
            kval(count)=k;
            count=count+1;
        end
        mean2ndk=sum(fdik.*kval);
        tempkmin=kmin_round;
        tempkmax=kmax;
        diff_observed_f(i,j)=mean2ndk-meank_f;
        kmax_round=kmax;
        while abs(mean2ndk-meank_f)>howclose % and i use the iterative search to find a kmax that minimizes the difference with fdik_orig
            clear simtillN 
            clear fdik
            clear kval
            kmax_round=round((tempkmax+tempkmin)/2);
            count=1;
            fdik=zeros(kmax_round-kmin_round,1);
            kval=zeros(kmax_round-kmin_round,1);
            kk=kmin_round:kmax_round;
            simtillN=sum((1./(kk.^gf)));
            for k=kmin_round:kmax_round
                fdik(count)=(1/k^(gf))/simtillN;
                kval(count)=k;
                count=count+1;
            end
            mean2ndk=sum(fdik.*kval);
            if tempkmax==tempkmin+1 || tempkmax==tempkmin; diff_observed_f(i,j)=mean2ndk-meank_f; flag_f(i,j)=1; break; end % this breaks if we reached an impasse.
            if mean2ndk-meank_f>0;  tempkmax=kmax_round; end
            if mean2ndk-meank_f<0; tempkmin=kmax_round; end
            diff_observed_f(i,j)=mean2ndk-meank_f;
        end

        fdik_orig=fdik;
        kval_orig_f=kval;
        kmin_orig_f=kmin_round;
        kmax_orig_f=kmax_round;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate baseline network g(k)
        kmin_round=kmin+1; % same as a above: i increase kmin of one to have a more condensed network
        count=1;
        gdik=zeros(kmax-kmin_round,1);
        kval=zeros(kmax-kmin_round,1);
        kk=kmin_round:kmax;
        simtillN=sum((1./(kk.^gg)));
        for k=kmin_round:1:kmax
            gdik(count)=(1/k^(gg))/simtillN;
            kval(count)=k;
            count=count+1;
        end
        mean2ndk=sum(gdik.*kval);
        tempkmin=kmin_round;
        tempkmax=kmax;
        diff_observed_g(i,j)=mean2ndk-meank_gg;
        kmax_round=kmax;
        while abs(mean2ndk-meank_gg)>howclose % and i use an iterative search to find a suitable kmax that minimzes the difference with fdik_orig
            clear simtillN 
            clear fdik
            clear kval
            kmax_round=round((tempkmax+tempkmin)/2);
            count=1;
            gdik=zeros(kmax_round-kmin,1);
            kval=zeros(kmax_round-kmin,1);
            kk=kmin_round:kmax_round;
            simtillN=sum((1./(kk.^gg)));
            for k=kmin_round:kmax_round
                gdik(count)=(1/k^(gg))/simtillN; 
                kval(count)=k;
                count=count+1;
            end
            mean2ndk=sum(gdik.*kval);
            if tempkmax==tempkmin+1 || tempkmax==tempkmin; diff_observed_g(i,j)=mean2ndk-meank_gg; flag_g(i,j)=1; break; end % this breaks if we reached an impasse.
            if mean2ndk-meank_gg>0;  tempkmax=kmax_round; end
            if mean2ndk-meank_gg<0; tempkmin=kmax_round; end
            diff_observed_g(i,j)=mean2ndk-meank_gg;
        end

        gdik_orig=gdik; % mi salvo k e fdik per il mps
        kval_orig_g=kval;
        kmin_orig_g=kmin_round;
        kmax_orig_g=kmax_round;

        
        % compute the betaa such that the consistency condition is
        % fulfilled.
        betaa(i,j)=meank_f/(meank_gg+meank_f);


        % compute effect on profit of mps on F
        [effect_f(i,j),klower_star_mps_f(i,j),L_star_mps_f(i,j),GAMMA_star_mps_f(i,j),phi_star_mps_f(i,j),b_star_mps_f(i,j),klower_star_orig_f(i,j),L_star_orig_f(i,j),GAMMA_star_orig_f(i,j),...
            phi_star_orig_f(i,j),b_star_orig_f(i,j),profit_mps_f(i,j),profit_orig_f(i,j)]=effect_of_mps_f(kmax_mps_f,kmin_mps_f,kmax_orig_f,kmin_orig_f,kval_mps_f,kval_orig_f,kval_orig_g,gdik_orig,fdik_orig,fdik_mps,...
            cost,betaa(i,j),p1,p2);

        % compute effect on profit of mps on G
        [effect_g(i,j),klower_star_mps_g(i,j),L_star_mps_g(i,j),GAMMA_star_mps_g(i,j),phi_star_mps_g(i,j),b_star_mps_g(i,j),klower_star_orig_g(i,j),L_star_orig_g(i,j),GAMMA_star_orig_g(i,j),...
            phi_star_orig_g(i,j),b_star_orig_g(i,j),profit_mps_g(i,j),profit_orig_g(i,j)]=effect_of_mps_g(kmax_mps_g,kmin_mps_g,kmax_orig_g,kmin_orig_g,kval_mps_g,kval_orig_g,kval_orig_f,gdik_orig,fdik_orig,gdik_mps,...
            cost,betaa(i,j),p1,p2);  

    end
end


% i print the results for the MPS in g(k)
figure(13)
imagesc(gammf,gammg,profit_mps_g-profit_orig_g)
set(gca,'YDir','normal')
xlabel('\gamma_f','FontSize',18);
ylabel('\gamma_g','FontSize',18);
colorbar
title('Effect of a MPS on g(k) on Difference in \pi*','Fontsize',20)
print -depsc SF_FdifferentG_pi_G.eps

% i print the results for the MPS in f(k)
figure(14)
imagesc(gammf,gammg,profit_mps_f-profit_orig_f)
set(gca,'YDir','normal')
xlabel('\gamma_f','FontSize',18);
ylabel('\gamma_g','FontSize',18);
colorbar
title('Effect of a MPS on f(k) on Difference in \pi*','Fontsize',20)
print -depsc SF_FdifferentG_pi_F.eps

    
% figure(1)
% imagesc(gammf,gammg,effect_f)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Sign of the effect of a MPS on f(k) on profits','Fontsize',20)
% print -depsc SF_FdifferentG_effect_F.eps
% 
% 
% figure(2)
% imagesc(gammf,gammg,klower_star_mps_f-klower_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on K*','Fontsize',20)
% print -depsc SF_FdifferentG_K_F.eps
% 
% 
% figure(3)
% imagesc(gammf,gammg,L_star_mps_f-L_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on L*','Fontsize',20)
% print -depsc SF_FdifferentG_L_F.eps
% 
% 
% figure(4)
% imagesc(gammf,gammg,GAMMA_star_mps_f-GAMMA_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on \Gamma*','Fontsize',20)
% print -depsc SF_FdifferentG_GAMMA_F.eps
% 
% 
% figure(5)
% imagesc(gammf,gammg,phi_star_mps_f-phi_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on \phi*','Fontsize',20)
% print -depsc SF_FdifferentG_phi_F.eps
% 
% figure(6)
% imagesc(gammf,gammg,b_star_mps_f-b_star_orig_f)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on f(k) on b*','Fontsize',20)
% print -depsc SF_FdifferentG_b_F.eps
% 
% 
% figure(7)
% imagesc(gammf,gammg,effect_g)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Sign of the effect of a MPS on g(k) on profits','Fontsize',20)
% print -depsc SF_FdifferentG_effect_G.eps
% 
% 
% figure(8)
% imagesc(gammf,gammg,klower_star_mps_g-klower_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on K*','Fontsize',20)
% print -depsc SF_FdifferentG_K_G.eps
% 
% figure(9)
% imagesc(gammf,gammg,L_star_mps_g-L_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on L*','Fontsize',20)
% print -depsc SF_FdifferentG_L_G.eps
% 
% 
% figure(10)
% imagesc(gammf,gammg,GAMMA_star_mps_g-GAMMA_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on \Gamma*','Fontsize',20)
% print -depsc SF_FdifferentG_GAMMA_G.eps
% 
% 
% figure(11)
% imagesc(gammf,gammg,phi_star_mps_g-phi_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on \phi*','Fontsize',20)
% print -depsc SF_FdifferentG_phi_G.eps
% 
% 
% figure(12)
% imagesc(gammf,gammg,b_star_mps_g-b_star_orig_g)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Effect of a MPS on g(k) on b*','Fontsize',20)
% print -depsc SF_FdifferentG_b_G.eps
% figure(15)
% imagesc(gammf,gammg,diff_observed_f)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Numerical Difference mean(ORIG_f)-mean(MPS_f)','Fontsize',20)
% print -depsc SF_FdifferentG_DiffObs_f.eps
% 
% figure(16)
% imagesc(gammf,gammg,diff_observed_g)
% set(gca,'YDir','normal')
% xlabel('\gamma_f','FontSize',18);
% ylabel('\gamma_g','FontSize',18);
% colorbar
% title('Numerical Difference mean(ORIG_g)-mean(MPS_g)','Fontsize',20)
% print -depsc SF_FdifferentG_DiffObs_g.eps
