%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NON-PARAMETRIC FRF:
% -------------------
% Descr.:   example of measurement data pre-processing
%           and non-parametric synchronized FRF estimation
% System:   Conventional motor-bench with flexible coupling
% Author:   Thomas Beauduin, KULeuven, PMA division, 2014
%           Wataru Ohnishi, The University of Tokyo, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('private/MultisineTypeA.mat');     % schoeder multisine experiment

%% STEP 1: TIME TREATMENT
% remove transient periods, offsets and trends
trans = 1;                      % number of transient periods
trend = 0;                      % period trend removal flag
r0 = (1:ms.nrofs);                 % time visualization range
rn = (ms.nrofs*2+1:(2+1)*ms.nrofs);   % data visualization range
input = [iq_adx];               % input data to motor bench
output = [theta_mx,-theta_my];  % output data of motor bench
[x,time] = pretreat(input,ms.nrofs,ms.harm.fs,trans,trend);
[y,time] = pretreat(output,ms.nrofs,ms.harm.fs,trans,trend);

figure
subplot(211), plot(time(r0),x(rn,:))
    title('input data'), legend('motor torque'), ylim([-20,20])
subplot(212), plot(time(r0),y(rn,:))
    title('output data'), legend('motor angle','load angle')
    xlabel('time [s]')
pause

%% STEP 2: NON-PARAMETRIC ESTIMATION
% fft data and vizualize in freq domain position data
flagTime = true;
Pest = time2frf_ml(x,y,ms,flagTime);

title('Motor-side');
bode_fdi({Pest(1,1)},[Pest.freq,Pest.UserData.FRFn(:,1)]);
legend('FRF','FRFn');

title('Load-side');
bode_fdi({Pest(2,1)},[Pest.freq,Pest.UserData.FRFn(:,2)]);
legend('FRF','FRFn');

% figure
% subplot(221), semilogx(Pest.freq,dbm(squeeze(Pest.resp(1,1,:))),Pest.freq,dbm(Pest.UserData.FRFn(:,1)),'r');
%     title('Motor-side'), ylabel('Magnitude [dB]');
% subplot(223), semilogx(Pest.freq,phs(squeeze(Pest.resp(1,1,:)),1))
%     xlabel('Frequency [Hz]'), ylabel('Phase [deg]');
% subplot(222), semilogx(Pest.freq,dbm(squeeze(Pest.resp(2,1,:))),Pest.freq,dbm(Pest.UserData.FRFn(:,2)),'r');
%     title('Load-side'), ylabel('Magnitude [dB]');
% subplot(224), semilogx(Pest.freq,phs(squeeze(Pest.resp(2,1,:)),1))
%     xlabel('Frequency [Hz]'), ylabel('Phase [deg]');
 
