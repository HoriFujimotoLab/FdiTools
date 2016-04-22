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

%% STEP 4: PARAMETRIC ESTIMATION
% deterministric/stochastic estimation with non-parametric noise model
n=4;                        % model order of denominator polynomial
mh=[2,0]; ml=[0,0];         % model orders of numerator polynomial
relvar=1e-10;               % relative variation of costfunction (stop)
iter=5e2;                   % maximum number of iterations (stop)
GN = 0;                     % Levenberg-Marquardt optimization
cORd = 'c';                 % continuous model identifaction
FRF_W = ones(size(FRFs));   % least squares weighting function
relax = 1;                  % relaxation factor for btls estimation

% Deterministic methods: non-linear least squares
[SYS.nls,SYS.wls] = nlsfdi(FRFs,freq,FRF_W,n,mh,ml,iter,relvar,GN,cORd,fs);
    FRF.nls = hfrf(SYS.nls,freq);
    FRF.wls = hfrf(SYS.wls,freq);
% Stochastic methods: maximum likelihood estimation
[SYS.ml,SYS.ls] = mlfdi(X,Y,freq,sX2,sY2,cXY,n,mh,ml,iter,relvar,GN,cORd,fs);
    FRF.ml = hfrf(SYS.ml,freq);
    FRF.ls = hfrf(SYS.ls,freq);

% Stochastic methods: bootstrapped total least squares
[SYS.btls,SYS.gtls] = btlsfdi(X,Y,freq,n,mh,ml,sY2,sX2,cXY,relax,iter,relvar,cORd,fs);
    FRF.btls = hfrf(SYS.btls,freq);
    FRF.gtls = hfrf(SYS.gtls,freq);

figure('Name','Deterministic Estimators')
subplot(221), semilogx(freq,[dbm(FRFs(:,1)),dbm(FRF.wls(:,1)),dbm(FRF.nls(:,1))])
hold on, semilogx(freq,dbm(FRF.ls(:,1)),'k--')    % starting values
hold on, semilogx(freq,dbm(FRFn(:,1)),'m--')      % uncertainty level
    ylabel('Amplitude [dB]') 
    legend('FRF','WLS','NLS','LS','FRFn')
    xlim([10,300])
subplot(223), semilogx(freq,[phs(FRFs(:,1)),phs(FRF.wls(:,1)),phs(FRF.nls(:,1))])
hold on, semilogx(freq,phs(FRF.ls(:,1)),'k--')    % starting values
    ylabel('Phase [deg]')
    xlabel('frequency [Hz]')
    xlim([10,300])
subplot(222), semilogx(freq,[dbm(FRFs(:,2)),dbm(FRF.wls(:,2)),dbm(FRF.nls(:,2))])
hold on, semilogx(freq,dbm(FRF.ls(:,2)),'k--')    % starting values
hold on, semilogx(freq,dbm(FRFn(:,2)),'m--')      % uncertainty level
    ylabel('Amplitude [dB]') 
    legend('FRF','WLS','NLS','LS','FRFn')
    xlim([10,300])
subplot(224), semilogx(freq,[phs(FRFs(:,2)),phs(FRF.wls(:,2)),phs(FRF.nls(:,2))])
hold on, semilogx(freq,phs(FRF.ls(:,2)),'k--')    % starting values
    ylabel('Phase [deg]')
    xlabel('frequency [Hz]')
    xlim([10,300])
    
figure('Name','Stochastic Estimators')
subplot(221), semilogx(freq,[dbm(FRFs(:,1)),dbm(FRF.ml(:,1)),dbm(FRF.btls(:,1))])
hold on, semilogx(freq,dbm(FRF.gtls(:,1)),'k--')  % starting values
hold on, semilogx(freq,dbm(sCR(:,1)),'m--')       % uncertainty lowerbound
    ylabel('Amplitude [dB]')
    legend('FRF','MLE','BTLS','GTLS','sCR')
    xlim([10,300])
subplot(223), semilogx(freq,[phs(FRFs(:,1)),phs(FRF.ml(:,1)),phs(FRF.btls(:,1))])
hold on, semilogx(freq,phs(FRF.gtls(:,1)),'k--')  % starting values
    ylabel('Phase [deg]')
    xlabel('frequency [Hz]')
    xlim([10,300])    
subplot(222), semilogx(freq,[dbm(FRFs(:,2)),dbm(FRF.ml(:,2)),dbm(FRF.btls(:,2))])
hold on, semilogx(freq,dbm(FRF.gtls(:,2)),'k--')  % starting values
hold on, semilogx(freq,dbm(sCR(:,2)),'m--')       % uncertainty lowerbound
    ylabel('Amplitude [dB]')
    legend('FRF','MLE','BTLS','GTLS','sCR')
    xlim([10,300])
subplot(224), semilogx(freq,[phs(FRFs(:,2)),phs(FRF.ml(:,2)),phs(FRF.btls(:,2))])
hold on, semilogx(freq,phs(FRF.gtls(:,2)),'k--')  % starting values
    ylabel('Phase [deg]')
    xlabel('frequency [Hz]')
    xlim([10,300])    
