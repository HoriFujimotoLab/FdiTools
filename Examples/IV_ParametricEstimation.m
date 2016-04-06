%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRIC ESTIMATION:
% ----------------------
% Descr.:   Example of parametric system model estimation
%           from frequency-domain data with known noise model.
% System:   Conventional motor-bench with flexible coupling.
% Author:   Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('MultisineTypeA.mat');     % schoeder multisine experiment

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
mh=[3,1]; ml=[0,0];         % model orders of numerator polynomial
relvar=1e-10;               % relative variation of costfunction (stop)
iter=5e2;                   % maximum number of iterations (stop)
GN = 0;                     % Levenberg-Marquardt optimization
cORd = 'c';                 % continuous model identifaction
FRF_W = ones(size(FRFs));   % least squares weighting function
relax = 1;                  % relaxation factor for btls estimation

% Deterministic methods: non-linear least squares
[Bn_nls,An_nls,Bn_wls,An_wls] = ...
nlsfdi(FRFs,freq,FRF_W,n,mh,ml,iter,relvar,GN,cORd,fs);
    SYS.nls(:,1) = tf(Bn_nls(1,:),An_nls); 
    FRF.nls(:,1) = squeeze(freqresp(SYS.nls(:,1),freq*2*pi));
    SYS.wls(:,1) = tf(Bn_wls(1,:),An_wls); 
    FRF.wls(:,1) = squeeze(freqresp(SYS.wls(:,1),freq*2*pi));
    SYS.nls(:,2) = tf(Bn_nls(2,:),An_nls); 
    FRF.nls(:,2) = squeeze(freqresp(SYS.nls(:,2),freq*2*pi));
    SYS.wls(:,2) = tf(Bn_wls(2,:),An_wls); 
    FRF.wls(:,2) = squeeze(freqresp(SYS.wls(:,2),freq*2*pi));
figure;
subplot(221),semilogx(freq,[dbm(FRFs(:,1)),dbm(FRF.wls(:,1)),dbm(FRF.nls(:,1))])
title('NLS MIMO')
subplot(223),semilogx(freq,[phs(FRFs(:,1)),phs(FRF.wls(:,1)),phs(FRF.nls(:,1))])
subplot(222),semilogx(freq,[dbm(FRFs(:,2)),dbm(FRF.wls(:,2)),dbm(FRF.nls(:,2))])
subplot(224),semilogx(freq,[phs(FRFs(:,2)),phs(FRF.wls(:,2)),phs(FRF.nls(:,2))])
