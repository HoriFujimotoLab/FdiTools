%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRIC ESTIMATION:
% ----------------------
% Descr.:   Example of parametric system model estimation
%           from frequency-domain data with known noise model.
% System:   Conventional motor-bench with flexible coupling.
% Author:   Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('private/MultisineTypeA.mat');     % schoeder multisine experiment

% Time treatment: remove transients/offsets/trends
trans = 1;                      % number of transient periods
trend = 0;                      % period trend removal flag
input = [iq_adx];               % input data to motor bench
output = [theta_mx,-theta_my];  % output data of motor bench
[x,time] = pretreat(input,nrofs,fs,trans,trend);
[y,time] = pretreat(output,nrofs,fs,trans,trend);

% Non-parametric estimation of frf matrix data
[X,Y,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,y,fs,fl,fh,df);

% deterministric/stochastic estimation with non-parametric noise model
n=4;                        % model order of denominator polynomial
mh=[2,0]; ml=[0,0];         % model orders of numerator polynomial
relvar=1e-10;               % relative variation of costfunction (stop)
iter=5e2;                   % maximum number of iterations (stop)
GN = 0;                     % Levenberg-Marquardt optimization
cORd = 'c';                 % continuous model identifaction
FRF_W = ones(size(FRFs));   % least squares weighting function
relax = 1;                  % relaxation factor for btls estimation

% Deterministic method: non-linear least squares
[Bn_nls,An_nls,Bn_wls,An_wls] = ...
nlsfdi(FRFs,freq,FRF_W,n,mh,ml,iter,relvar,GN,cORd,fs);
    SYS.nls(:,1) = tf(Bn_nls(1,:),An_nls); 
    SYS.wls(:,1) = tf(Bn_wls(1,:),An_wls); 
    SYS.nls(:,2) = tf(Bn_nls(2,:),An_nls);
    SYS.wls(:,2) = tf(Bn_wls(2,:),An_wls); 
    
% Stochastic method: maximum likelihood estimation
[Bn_ml,An_ml,Bn_ls,An_ls] = ...
mlfdi(X,Y,freq,sX2,sY2,cXY,n,mh,ml,iter,relvar,GN,cORd,fs);
    SYS.ml(:,1) = tf(Bn_ml(1,:),An_ml);
    SYS.ls(:,1) = tf(Bn_ls(1,:),An_ls);
    SYS.ml(:,2) = tf(Bn_ml(2,:),An_ml);
    SYS.ls(:,2) = tf(Bn_ls(2,:),An_ls);

% Stochastic method: bootstrapped total least squares
[Bn_btls,An_btls,Bn_gtls,An_gtls] = ...
btlsfdi(X,Y,freq,n,mh,ml,sY2,sX2,cXY,relax,iter,relvar,cORd,fs);
    SYS.btls(:,1) = tf(Bn_btls(1,:),An_btls); 
    SYS.gtls(:,1) = tf(Bn_gtls(1,:),An_gtls);
    SYS.btls(:,2) = tf(Bn_btls(2,:),An_btls); 
    SYS.gtls(:,2) = tf(Bn_gtls(2,:),An_gtls);

% PART 5: SELECTION & VALIDATION
% Selection & validation of estimator type and model set
nrofp = length(x)/fs*df;        % number of measured periods

%% TEST 1: RESIDUALS
% whiteness test with fraction above the x%-percentile confidence bound
[lags,corr,cb50,frac50,tag,cb95,frac95] = residtest(x,y,freq,FRFs,SYS,sCR,fs);
figure 
plot(lags,corr(:,:,1),'.',lags,cb95,'k',lags,cb50,'k--')
    title(strcat(tag(:,1),' -- p_{>95%}:',num2str(frac95(:,1)),...
          ' _ - _ p_{>50%}:',num2str(frac50(:,1))));
    legend([tag(:,1);'frac_{95%}';'frac_{50%}'])
    xlabel('Lag number'), ylabel('Amplitude [dB]'), ylim([0,1.5])

figure
plot(lags,corr(:,:,2),'.',lags,cb95,'k',lags,cb50,'k--')
    title(strcat(tag(:,1),' -- p_{>95%}:',num2str(frac95(:,2)),...
          ' _ - _ p_{>50%}:',num2str(frac50(:,2))));
    legend([tag(:,1);'frac_{95%}';'frac_{50%}'])
    xlabel('Lag number'), ylabel('Amplitude [dB]'),ylim([0,1.5])
pause

%% TEST 2: COST FUNCTION
% residual estimation cost (maximum likelihood function)
[cost,intv,tag] = costtest(X,Y,freq,sX2,sY2,cXY,SYS,relax,nrofp);
figure
bar(cost,0.8), hold on
line([0,length(cost)+1],[intv(2),intv(2)],'Color','k','LineStyle','--')
    set(gca,'xticklabel',tag)
    title('Estimator Selection: Residual Cost')
    legend('H_{11}','H_{12}','noise'), ylim([0,intv(2)*5])
pause

%% TEST 3: CHI-SQUARES
% Chi^2 test of absolute modeling error
[confid,var,tag] = chi2test(X,Y,freq,FRFs,sCR,SYS);
figure
subplot(211),loglog(freq,confid(:,:,1)), hold on, loglog(freq,var(:,1),'k')
    legend(tag(:,1),'crlb')
subplot(212),loglog(freq,confid(:,:,2)), hold on, loglog(freq,var(:,2),'k')
    legend(tag(:,2),'crlb') 
    
