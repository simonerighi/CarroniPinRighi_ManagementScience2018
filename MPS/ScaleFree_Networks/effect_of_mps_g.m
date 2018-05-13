function [effect,klower_star_mps,L_star_mps,GAMMA_star_mps,phi_star_mps,b_star_mps,klower_star_orig,L_star_orig,GAMMA_star_orig,...
    phi_star_orig,b_star_orig,pistar_mps,pistar_orig]=effect_of_mps_g(kmax_mps,kmin_mps,kmax_orig,kmin_orig,kval_mps,kval_orig_g,kval_orig_f,gdik_orig,fdik_orig,gdik_mps,...
    cost,betaa,p1,p2)
% this function computes the maximum payoff given the input k f(k) and
% g(k) both in baseline and in the case of a Mean Preserving Spread
% (Second order stochastic dominant shift). The function passes as output
% the optimal values of all key variables of the model (chiefly the profit
% pi) as well as the comparison between the two profits in the variable
% "effect".
% this function works for MPS performed on the distribution g(k).


% for all k compute the value of profit and other key variables in the MPS network 

L=zeros(kmax_mps-kmin_mps+1,1);
GAMMA=zeros(kmax_mps-kmin_mps+1,1);
phi=zeros(kmax_mps-kmin_mps+1,1);
b=zeros(kmax_mps-kmin_mps+1,1);
pi=zeros(kmax_mps-kmin_mps+1,1);


count=1;
for klower=kmin_orig:1:kmax_orig
    L(count)=(1-p1)*sum(fdik_orig(find(kval_orig_f==klower):end)); 
    GAMMA(count)=1-sum(gdik_mps.*((1-L(count)).^kval_mps));
    phi(count)=(1-p2)*sum(gdik_mps.*((1-(1-L(count)).^kval_mps)./(kval_mps.*L(count))));
    b(count)=cost/(phi(count)*klower);
    pi(count)=(1-betaa)*p1*(1-p1)+betaa*(p2-b(count))*(1-p2)*GAMMA(count);
    count=count+1;
end

% find the maximum profit in the MPS
[pistar_mps, IX]=max(pi); % maximum profit
klower_star_mps=IX+kmin_mps-1; % klowerbar_star
L_star_mps=L(IX);
GAMMA_star_mps=GAMMA(IX);
phi_star_mps=phi(IX);
b_star_mps=b(IX);


% for all k compute the value of profit and other key variables in the baseline network 

clear L
clear GAMMA
clear phi
clear b
clear pi

 L=zeros(kmax_orig-kmin_orig+1,1);
 GAMMA=zeros(kmax_orig-kmin_orig+1,1);
 phi=zeros(kmax_orig-kmin_orig+1,1);
 b=zeros(kmax_orig-kmin_orig+1,1);
 pi=zeros(kmax_orig-kmin_orig+1,1);

count=1;
for klower=kmin_orig:1:kmax_orig
    L(count)=(1-p1)*sum(fdik_orig(find(kval_orig_f==klower):end)); % probability that an informed buyer passes the info.
    GAMMA(count)=1-sum(gdik_orig.*((1-L(count)).^kval_orig_g));  
    phi(count)=(1-p2)*sum(gdik_orig.*((1-(1-L(count)).^kval_orig_g)./(kval_orig_g.*L(count))));
    b(count)=cost/(phi(count)*klower);
    pi(count)=(1-betaa)*p1*(1-p1)+betaa*(p2-b(count))*(1-p2)*GAMMA(count);
    count=count +1;
end

% find the maximum profit in the baseline network

[pistar_orig, IX]=max(pi); % maximum profit
klower_star_orig=IX+kmin_orig-1; % klowerbar_star
L_star_orig=L(IX);
GAMMA_star_orig=GAMMA(IX);
phi_star_orig=phi(IX);
b_star_orig=b(IX);



%the MPS shift...
if pistar_orig>pistar_mps 
    effect=-1; % ...decreases...
end

 if pistar_orig<pistar_mps
     effect=1; %...increases...
 end
 
 if pistar_orig==pistar_mps
     effect=0; %...leaves equal...
 end
 %..profits