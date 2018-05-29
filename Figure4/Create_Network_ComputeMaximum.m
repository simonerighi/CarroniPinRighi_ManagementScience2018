function [maxpi_orig,maxL_orig,maxGAMMA_orig,maxphi_orig,klowerbest,klowerbestrelative,klowerbest_only_public]=Create_Network_ComputeMaximum(n,gammag1,betaa)
% Generates results for f being Scale Free and g being random.
% Carroni Pin Righi (2018), submitted for pubblication in Management
% Sciences


p1=1/2;
p2=1/2;
cost=0.3;
kmin=1;
kmax=100;
warning off
gammg=gammag1;
howclose=1.0000e-13;

gg=gammg;
%%%%% %%%%% %%%%% %%%%% 
%%%%% OUTDEGREE f (SF)
%%%%% %%%%% %%%%% %%%%% 

%network pdf creation
kval=[1:1:n-2];
fdik=zeros(kmax-kmin+1,1);
kval=zeros(kmax-kmin+1,1);
kk=kmin:kmax;
simtillN=sum((1./(kk.^gg)));
count=1;
for k=kmin:1:kmax
    fdik(count)=(1/k^(gg))/simtillN;
    kval(count)=k;
    count=count+1;
end
mean_fsearch=sum(fdik.*kval);

fdik_orig=fdik; % mi salvo k e fdik per il mps
kval_orig_f=kval;
kmin_orig_f=kmin;
kmax_orig_f=kmax;



fcorr=0.5;

gg=((mean_fsearch)/n)*((betaa)/(1-betaa));  % + (1.1/n);



%%%%% %%%%% %%%%% %%%%% 
%%%%% INDEGREE g (RANDOM)
%%%%% %%%%% %%%%% %%%%% 


kval=[1:1:n-2];
gdik=binopdf(kval,n-2,gg);
meank_g=sum(gdik.*kval);

mingamma=0.01;
maxgamma=1;

diff_observed=mean_fsearch-meank_g;
countc=1;
while abs(diff_observed)>howclose % ricerca iterativa di un valore di kmax che minimizzi la differenza con fdik_orig
    clear simtillN
    clear gdik
    
    if diff_observed<0;  mingamma=gg; end
    if diff_observed>0;  maxgamma=gg; end    
    gg=mingamma+((maxgamma-mingamma)/2);
    
    % compute new network with gg slope
    gdik=binopdf(kval,n-2,gg);
    
    meank_g=sum(gdik.*kval);
    diff_observed=mean_fsearch-meank_g;
    countc=countc+1;
    if countc>=200 
        set_to_nan; 
        return % if it can't find the the right gg in 200 rounds the software will return error.
    end
end





%I eliminate extreme values for which the pdf assumes value 0 (if
%necessary)
dummygdik=gdik;
if find(dummygdik==0,1)==1; kmin_orig_g=find(dummygdik>0,1); dummygdik(1:kmin_orig_g-1)=dummygdik(1:kmin_orig_g-1)+1; else; kmin_orig_g=1; end;
if dummygdik(end)==0; kmax_orig_g=find(dummygdik==0,1)-1; else; kmax_orig_g=n-2; end; %length(dummyfdik)
clear dummyfdik
gdik_orig=gdik(kmin_orig_g:kmax_orig_g); % mi salvo k e fdik per il mps
kval_orig_g=kmin_orig_g:1:kmax_orig_g;
clear simtillN 
clear gdik
clear kval




%%%%%%%%%%%%%%%%%%%%
%%%% Compute variable values for all k
%%%%%%%%%%%%%%%%%%%%

Lorig=zeros(kmax_orig_f-kmin_orig_f+1,1);
GAMMA_orig=zeros(kmax_orig_f-kmin_orig_f+1,1);
phi_orig=zeros(kmax_orig_f-kmin_orig_f+1,1);
B=zeros(kmax_orig_f-kmin_orig_f+1,1);
pi_orig=zeros(kmax_orig_f-kmin_orig_f+1,1);


count=1;
p2=1/2;

for klower=kmin_orig_f:1:kmax_orig_f
    Lorig(count)=sum(fdik_orig(find(kval_orig_f==klower):end)).*(1+sum(fdik_orig(1:find(kval_orig_f==klower)))); % probability that an informed buyer passes the info.
    GAMMA_orig(count)=1-sum(gdik_orig.*((1-Lorig(count)).^kval_orig_g));
    phi_orig(count)=(1-p2)*sum(gdik_orig.*((1-(1-Lorig(count)).^kval_orig_g)./(kval_orig_g.*Lorig(count))));
    if klower > kmin_orig_f % this solves the problem of numerical precision
        if phi_orig(count-1)>0.49999999
            phi_orig(count)=0.5;
        end
    end
    B(count)=(cost./phi_orig(count)).*sum(fdik_orig(find(kval_orig_f==klower):end)).*((1+klower)/(2*klower) + sum(fdik_orig(1:find(kval_orig_f==klower))));
    pi_orig(count)=(p2-B(count)).*GAMMA_orig(count);
    count=count+1;
end

% find maximum value and associated position.
[maxpi_orig,locmax_orig]=max(pi_orig);
maxL_orig=Lorig(locmax_orig);
maxGAMMA_orig=GAMMA_orig(locmax_orig);
maxphi_orig=phi_orig(locmax_orig);
klowerbest=kmin_orig_f+locmax_orig;
klowerbestrelative=(kmin_orig_f+locmax_orig-1)/kmax_orig_f;


% This gives the optimal profits for the case in which, on the same
% network, only public WOM is considered

Lorig=zeros(kmax_orig_f-kmin_orig_f+1,1);
GAMMA_orig=zeros(kmax_orig_f-kmin_orig_f+1,1);
phi_orig=zeros(kmax_orig_f-kmin_orig_f+1,1);
b_orig=zeros(kmax_orig_f-kmin_orig_f+1,1);
pi_orig=zeros(kmax_orig_f-kmin_orig_f+1,1);


count=1;
for klower=kmin_orig_f:1:kmax_orig_f
    Lorig(count)=(1-p1)*sum(fdik_orig(find(kval_orig_f==klower):end)); % probability that an informed buyer passes the info.
    GAMMA_orig(count)=1-sum(gdik_orig.*((1-Lorig(count)).^kval_orig_g));
    phi_orig(count)=(1-p2)*sum(gdik_orig.*((1-(1-Lorig(count)).^kval_orig_g)./(kval_orig_g.*Lorig(count))));
    if klower > kmin_orig_f % this solves the problem of numerical precision
        if phi_orig(count-1)>0.49999999
            phi_orig(count)=0.5;
        end
    end
    b_orig(count)=cost/(phi_orig(count)*klower);
    
    pi_orig(count)=(1-betaa(1))*p1*(1-p1)+betaa(1)*(p2-b_orig(count))*(1-p2)*GAMMA_orig(count);
    count=count+1;
end

[maxpi_orig_2,locmax_orig_2]=max(pi_orig);
maxL_orig_2=Lorig(locmax_orig_2);
maxGAMMA_orig_2=GAMMA_orig(locmax_orig_2);
maxphi_orig_2=phi_orig(locmax_orig_2);
maxb_orig_2=b_orig(locmax_orig_2);
klowerbest_only_public=kmin_orig_f+locmax_orig_2;