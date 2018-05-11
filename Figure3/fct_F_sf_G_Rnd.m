function [maxpi_orig,locmax_orig,maxL_orig,maxGAMMA_orig,maxphi_orig,maxb_orig,maxpi_FOSD,maxL_FOSD,maxGAMMA_FOSD,...
    maxphi_FOSD,maxb_FOSD,is_FOSD_Larger,Welfare_orig,Welfare_FOSD,diff_avg_orig,diff_avg_FOSD]=fct_F_sf_G_Rnd(n,gammag1,gammagfosd,betaa)
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
gammg=[gammag1 gammagfosd];
howclose=1.0000e-13;

%%%%% %%%%% %%%%% %%%%% 
%%%%% ORIGINALE - INDEGREE g (RANDOM)
%%%%% %%%%% %%%%% %%%%% 

gg=gammg(1); % to have only one variable

kval=[1:1:n-2];
gdik=binopdf(kval,n-2,gg);
meank_g=sum(gdik.*kval);

dummygdik=gdik;
if find(dummygdik==0,1)==1; kmin_orig_g=find(dummygdik>0,1); dummygdik(1:kmin_orig_g-1)=dummygdik(1:kmin_orig_g-1)+1; else; kmin_orig_g=1; end;
if dummygdik(end)==0; kmax_orig_g=find(dummygdik==0,1)-1; else; kmax_orig_g=n-2; end; %length(dummyfdik)
clear dummyfdik
gdik_orig=gdik(kmin_orig_g:kmax_orig_g); % mi salvo k e fdik per il mps
kval_orig_g=kmin_orig_g:1:kmax_orig_g;
clear simtillN 
clear gdik
clear kval


%%%%% %%%%% %%%%% %%%%% 
%%%%% ORIGINALE - OUTDEGREE f (SF)
%%%%% %%%%% %%%%% %%%%% 
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

gf=gg; % SF exponent for g is initalized as equal to the other;
maxgamma=10;
mingamma=1.001;

fmean_theo=(meank_g*(betaa)/(1-betaa)); %  value that the average g SHOULD assume to equalize the condition
diff_observed=mean_fsearch-fmean_theo;
countc=1;
while abs(diff_observed)>howclose % ricerca iterativa di un valore di kmax che minimizzi la differenza con fdik_orig
    clear simtillN
    clear gdik
    
    if diff_observed<0;  maxgamma=gf; end
    if diff_observed>0;  mingamma=gf; end    
    gf=mingamma+((maxgamma-mingamma)/2);
    
    % compute new network with gg slope
    count=1;
    fdik=zeros(kmax-kmin+1,1);
    kk=kmin:kmax;
    simtillN=sum((1./(kk.^gf)));
    for k=kmin:kmax
        fdik(count)=(1/k^(gf))/simtillN;
        kval(count)=k;
        count=count+1;
    end
    mean_fsearch=sum(fdik.*kval);
    countc=countc+1;
    diff_observed=mean_fsearch-fmean_theo; 
    if countc>=200; 
        set_to_nan; 
        return
    end;
end

fdik_orig=fdik; % mi salvo k e fdik per il mps
kval_orig_f=kval;
kmin_orig_f=kmin;
kmax_orig_f=kmax;


diff_avg_orig=((1-betaa)*mean_fsearch)-(betaa*meank_g);   


   
%%%%% %%%%% %%%%% %%%%% 
%%%%% %%%%% %%%%% %%%%% 
%%%%% %%%%% %%%%% %%%%% 

% GENERO FOSD g

gg=gammg(1); % to have only one variable

kval=[1:1:n-2];
gdik=binopdf(kval,n-2,gg);
meank_g2=sum(gdik.*kval);
mean_g_FOSD_mtx(1)=meank_g2;

dummygdik=gdik;
if find(dummygdik==0,1)==1; kmin_FOSD_g=find(dummygdik>0,1); dummygdik(1:kmin_FOSD_g-1)=dummygdik(1:kmin_FOSD_g-1)+1; else; kmin_FOSD_g=1; end;
if dummygdik(end)==0; kmax_FOSD_g=find(dummygdik==0,1)-1; else; kmax_FOSD_g=n-2; end; %length(dummyfdik)
clear dummyfdik
gdik_FOSD=gdik(kmin_FOSD_g:kmax_FOSD_g); % mi salvo k e fdik per il mps
kval_FOSD_g=kmin_FOSD_g:1:kmax_FOSD_g;
clear simtillN 
clear gdik
clear kval


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FOSD - OUTDEGREE SF... con beta fissato e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERO FOSD OUTEGREE f
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

gf=gg; % SF exponent for g is initalized as equal to the other;
maxgamma=10;
mingamma=1.01;

fmean_theo=meank_g2*(betaa/(1-betaa)); %  value that the average f SHOULD assume to equalize the condition
diff_observed=mean_fsearch-fmean_theo;
countc=1;
while abs(diff_observed)>howclose % ricerca iterativa di un valore di kmax che minimizzi la differenza con fdik_orig
    
    clear simtillN 
    clear fdik
    
    if mean_fsearch-fmean_theo<0;  maxgamma=gf; end
    if mean_fsearch-fmean_theo>0;  mingamma=gf; end    
    gf=mingamma+((maxgamma-mingamma)/2);
    
    count=1;
    fdik=zeros(kmax-kmin+1,1);
    simtillN=sum((1./(kk.^gf)));
    for k=kmin:kmax
        fdik(count)=(1/k^(gf))/simtillN;
        kval(count)=k;
        count=count+1;
    end
    mean_fsearch=sum(fdik.*kval);

    diff_observed=mean_fsearch-fmean_theo;
    countc=countc+1;
    if countc>=200; 
        set_to_nan; 
        return
    end;
end

fdik_FOSD=fdik; % mi salvo k e fdik per il mps
kval_FOSD_f=kval;
kmin_FOSD_f=kmin;
kmax_FOSD_f=kmax;

clear simtillN 
clear fdik
clear kval


diff_avg_FOSD=((1-betaa)*mean_fsearch)-(betaa*meank_g2);





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

if maxpi_FOSD>maxpi_orig
    is_FOSD_Larger=1;
else
    is_FOSD_Larger=0;
end

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


     
