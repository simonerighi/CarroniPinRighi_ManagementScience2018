function [maxpi_orig,locmax_orig,maxL_orig,maxGAMMA_orig,maxphi_orig,maxb_orig,maxpi_FOSD,maxL_FOSD,maxGAMMA_FOSD,...
    maxphi_FOSD,maxb_FOSD,is_FOSD_Larger,Welfare_orig,Welfare_FOSD,diff_avg_orig,diff_avg_FOSD]=fct_F_Rnd_G_Rnd(n,lambdag,lambdagfosd,betaa)
% genera risultati per f sf e g random.

%n=10000;
%gammf=[2.3 2.20];
%betaa=0.5;

p1=1/2;
p2=1/2;
cost=0.06;
kmin=1;
kmax=n-2;
delta_gradi=30; % this could be chosen relative to var(fdik*k);
warning off
lambdag=[lambdag lambdagfosd];

%%%%% %%%%% %%%%% %%%%% 
%%%%% ORIGINALE - INDEGREE g (RND)
%%%%% %%%%% %%%%% %%%%% 

gg=lambdag(1); % to have only one variable

kval=[1:1:n-2];
gdik=binopdf(kval,n-2,gg);

dummygdik=gdik;
if find(dummygdik==0,1)==1; kmin_orig_g=find(dummygdik>0,1); dummygdik(1:kmin_orig_g-1)=dummygdik(1:kmin_orig_g-1)+1; else; kmin_orig_g=1; end;
if dummygdik(end)==0; kmax_orig_g=find(dummygdik==0,1)-1; else; kmax_orig_g=n-2; end; %length(dummyfdik)
clear dummygdik
gdik_orig=gdik(kmin_orig_g:kmax_orig_g); % mi salvo k e fdik per il mps
kval_orig_g=kmin_orig_g:1:kmax_orig_g;
clear simtillN 
clear gdik
clear kval

meank_g=sum(gdik_orig.*kval_orig_g);

fcorr=0.5; %((10.1)/n)+gg*(betaa); %*((betaa)/(1-betaa));

gf=((meank_g+fcorr)/n)*((betaa)/(1-betaa));  % + (1.1/n);

%%%%% %%%%% %%%%% %%%%% 
%%%%% ORIGINALE - OUTDEGREE f (RANDOM)
% s.t. lambda=gg=(meank_f/n)*(betaa/(1-betaa));
%%%%% %%%%% %%%%% %%%%% 

kval=[1:1:n-2];
fdik=binopdf(kval,n-2,gf);

dummyfdik=fdik;
if find(dummyfdik==0,1)==1; kmin_orig_f=find(dummyfdik>0,1); dummyfdik(1:kmin_orig_f-1)=dummyfdik(1:kmin_orig_f-1)+1; else; kmin_orig_f=1; end;
if dummyfdik(end)==0; kmax_orig_f=find(dummyfdik==0,1)-1; else; kmax_orig_f=n-2; end; %length(dummyfdik)
clear dummyfdik
fdik_orig=fdik(kmin_orig_f:kmax_orig_f); % mi salvo k e fdik per il mps
kval_orig_f=kmin_orig_f:1:kmax_orig_f;
clear simtillN 
clear fdik
clear kval

meank_f=sum(fdik_orig.*kval_orig_f);


diff_avg_orig=((1-betaa)*(meank_f))-(betaa*(meank_g));   

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FOSD - INDEGREE Rand... con beta fissato e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gg=lambdag(2);

kval=[1:1:n-2];
gdik=binopdf(kval,n-2,gg);

dummygdik=gdik;
if find(dummygdik==0,1)==1; kmin_FOSD_g=find(dummygdik>0,1); dummygdik(1:kmin_FOSD_g-1)=dummygdik(1:kmin_FOSD_g-1)+1; else; kmin_FOSD_g=1; end;
if dummygdik(end)==0; kmax_FOSD_g=find(dummygdik==0,1)-1; else; kmax_FOSD_g=n-2; end; %length(dummyfdik)
clear dummygdik
gdik_FOSD=gdik(kmin_FOSD_g:kmax_FOSD_g); % mi salvo k e fdik per il mps
kval_FOSD_g=kmin_FOSD_g:1:kmax_FOSD_g;
clear simtillN 
clear gdik
clear kval

meank_g2=sum(gdik_FOSD.*kval_FOSD_g);

fcorr=0.5; %((10.1)/n)+gg*(betaa); %((1.1)/n)*();

% computo lambda da applicare
gf=((meank_g2+fcorr)/n)*((betaa(1))/(1-betaa(1))); % + (1.1/n);

%%%%%%%%

% GENERO FOSD f
kval=[1:1:n-2];
fdik=binopdf(kval,n-2,gf);
meank_f2=sum(fdik.*kval);

dummyfdik=fdik;
if find(dummyfdik==0,1)==1; kmin_FOSD_f=find(dummyfdik>0,1); dummyfdik(1:kmin_FOSD_f-1)=dummyfdik(1:kmin_FOSD_f-1)+1; else; kmin_FOSD_f=1; end;
if dummyfdik(end)==0; kmax_FOSD_f=find(dummyfdik==0,1)-1; else; kmax_FOSD_f=n-2; end; %length(dummyfdik)
clear dummyfdik
fdik_FOSD=fdik(kmin_FOSD_f:kmax_FOSD_f); % mi salvo k e fdik per il mps
kval_FOSD_f=kmin_FOSD_f:1:kmax_FOSD_f;
clear simtillN 
clear fdik
clear kval

diff_avg_FOSD=((1-betaa)*(meank_f2))-(betaa*(meank_g2));


%%%%%%%%%%%%%%%%%%%%
%%%% Graph for orig
%%%%%%%%%%%%%%%%%%%%

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
    if klower > kmin_orig_f; % this solves the problem of numerical precision
        if phi_orig(count-1)>0.49999999;
            phi_orig(count)=0.5;
        end
    end
    b_orig(count)=cost/(phi_orig(count)*klower);
    
    pi_orig(count)=(1-betaa(1))*p1*(1-p1)+betaa(1)*(p2-b_orig(count))*(1-p2)*GAMMA_orig(count);
    count=count+1;
end


%%%%%%%%%%%%%%%%%%%%
%%%% Graph for FOSD
%%%%%%%%%%%%%%%%%%%%

L=zeros(kmax_FOSD_f-kmin_FOSD_f+1,1);
GAMMA=zeros(kmax_FOSD_f-kmin_FOSD_f+1,1);
phi=zeros(kmax_FOSD_f-kmin_FOSD_f+1,1);
b=zeros(kmax_FOSD_f-kmin_FOSD_f+1,1);
pi_FOSD=zeros(kmax_FOSD_f-kmin_FOSD_f+1,1);


count=1;
for klower=kmin_FOSD_f:1:kmax_FOSD_f
    L(count)=(1-p1)*sum(fdik_FOSD(find(kval_FOSD_f==klower):end)); % probability that an informed buyer passes the info.
    GAMMA(count)=1-sum(gdik_FOSD.*((1-L(count)).^kval_FOSD_g));
    phi(count)=(1-p2)*sum(gdik_FOSD.*((1-(1-L(count)).^kval_FOSD_g)./(kval_FOSD_g.*L(count))));
    if klower > kmin_FOSD_f; % this solves the problem of numerical precision
        if phi(count-1)>0.49999999;
            phi(count)=0.5;
        end
    end
    b(count)=cost/(phi(count)*klower);
    pi_FOSD(count)=(1-betaa(1))*p1*(1-p1)+betaa(1)*(p2-b(count))*(1-p2)*GAMMA(count);
    count=count+1;
end

[maxpi_orig,locmax_orig]=max(pi_orig);
maxL_orig=Lorig(locmax_orig);
maxGAMMA_orig=GAMMA_orig(locmax_orig);
maxphi_orig=phi_orig(locmax_orig);
maxb_orig=b_orig(locmax_orig);

[maxpi_FOSD,locmax_FOSD]=max(pi_FOSD);
maxL_FOSD=L(locmax_FOSD);
maxGAMMA_FOSD=GAMMA(locmax_FOSD);
maxphi_FOSD=phi(locmax_FOSD);
maxb_FOSD=b(locmax_FOSD);

is_FOSD_Larger=maxpi_FOSD>maxpi_orig;

Welfare_orig=(3/8)*((1-betaa(1))+betaa(1)*maxGAMMA_orig);
Welfare_FOSD=(3/8)*((1-betaa(1))+betaa(1)*maxGAMMA_FOSD);





% 
% figure(20)
% semilogx(kmin_orig_f:1:kmax_orig_f,pi_orig,'b');
% hold on
% semilogx(kmin_FOSD_f:1:kmax_FOSD_f,pi_FOSD,'r');
% title('Indegree is Random, Outdegree SF')
% xlabel('K_{lowerbar}','FontSize',18)
% ylabel('Profits','FontSize',18)
% legend('\pi Original','\pi FOSD',0)
% saveas(20,'IndegreeRandomOutdegreeSF_FOSD','fig');


     
