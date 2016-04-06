%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NON-PARAMETRIC FRF:
% -------------------
% Descr.:   example of measurement data pre-processing
%           and non-parametric synchronized FRF estimation
% System:   Conventional motor-bench with flexible coupling
% Author:   Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('MultisineTypeA.mat');     % schoeder multisine experiment

%% STEP 1: TIME TREATMENT
% remove transient periods, offsets and trends
trans = 1;                      % number of transient periods
trend = 0;                      % period trend removal flag
r0 = (1:nrofs);                 % time visualization range
rn = (nrofs*5+1:(5+1)*nrofs);   % data visualization range
input = [iq_adx];               % input data to motor bench
output = [theta_mx,-theta_my];  % output data of motor bench
[x,time] = pretreat(input,nrofs,fs,trans,trend);
[y,time] = pretreat(output,nrofs,fs,trans,trend);

figure
subplot(211), plot(time(r0),x(rn,:))
    title('input data'), legend('motor torque'), ylim([-20,20])
subplot(212), plot(time(r0),y(rn,:))
    title('output data'), legend('motor angle','load angle')
    xlabel('time [s]')
pause

%% STEP 2: NON-PARAMETRIC ESTIMATION
% fft data and vizualize in freq domain position data
[X,Y,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,y,fs,fl,fh,df);

figure
subplot(221), semilogx(freq,dbm(FRFs(:,1)),freq,dbm(FRFn(:,1)),'r');
    title('Motor-side'), ylabel('Magnitude [dB]'), xlim([fl,fh])
subplot(223), semilogx(freq,phs(FRFs(:,1),1))
    xlabel('Frequency [Hz]'), ylabel('Phase [deg]'), xlim([fl,fh])
subplot(222), semilogx(freq,dbm(FRFs(:,2)),freq,dbm(FRFn(:,2)),'r');
    title('Load-side'), ylabel('Magnitude [dB]'), xlim([fl,fh])
subplot(224), semilogx(freq,phs(FRFs(:,2)))
    xlabel('Frequency [Hz]'), ylabel('Phase [deg]'), xlim([fl,fh])
 