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
input = [iq_adx];               % input data to motor bench
output = [theta_mx,-theta_my];  % output data of motor bench
[x,time] = pretreat(input,ms.nrofs,ms.harm.fs,trans,trend);
[y,time] = pretreat(output,ms.nrofs,ms.harm.fs,trans,trend);

nrofp = length(input)/ms.nrofs - trans; % number of period after transient removal

figure;
subplot(311)
for k = 1:nrofp
    h = plot(time(1:ms.nrofs),x(ms.nrofs*(k-1)+1:k*ms.nrofs)); hold on;
    h.DisplayName = sprintf('period%d',k+trans);
end
ylabel('input current [A]');
title(sprintf('%d periods (%d transient removal)',nrofp,trans));
% legend;
subplot(312)
for k = 1:nrofp
    h = plot(time(1:ms.nrofs),y(ms.nrofs*(k-1)+1:k*ms.nrofs,1)); hold on;
    h.DisplayName = sprintf('period%d',k+trans);
end
ylabel('motor angle [rad]');
subplot(313)
for k = 1:nrofp
    h = plot(time(1:ms.nrofs),y(ms.nrofs*(k-1)+1:k*ms.nrofs,2)); hold on;
    h.DisplayName = sprintf('period%d',k+trans);
end
ylabel('load angle [rad]');
xlabel('time [s]');

pause

%% STEP 2: NON-PARAMETRIC ESTIMATION
% fft data and vizualize in freq domain position data
flagTime = true;
Pest = time2frf_ml(x,y,ms,flagTime);

title('Motor-side');
bode_fdi({Pest(1,1)},[Pest.freq,Pest.UserData.sGhat(:,1)]);
legend('FRF','sGhat');

title('Load-side');
bode_fdi({Pest(2,1)},[Pest.freq,Pest.UserData.sGhat(:,2)]);
legend('FRF','sGhat');

% figure
% subplot(221), semilogx(Pest.freq,dbm(squeeze(Pest.resp(1,1,:))),Pest.freq,dbm(Pest.UserData.FRFn(:,1)),'r');
%     title('Motor-side'), ylabel('Magnitude [dB]');
% subplot(223), semilogx(Pest.freq,phs(squeeze(Pest.resp(1,1,:)),1))
%     xlabel('Frequency [Hz]'), ylabel('Phase [deg]');
% subplot(222), semilogx(Pest.freq,dbm(squeeze(Pest.resp(2,1,:))),Pest.freq,dbm(Pest.UserData.FRFn(:,2)),'r');
%     title('Load-side'), ylabel('Magnitude [dB]');
% subplot(224), semilogx(Pest.freq,phs(squeeze(Pest.resp(2,1,:)),1))
%     xlabel('Frequency [Hz]'), ylabel('Phase [deg]');
 
