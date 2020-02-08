%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NL DISTORTIONS DETECTION:
% -------------------------
% Descr.:   example of non-linear distortions detection
%           using odd-odd random phase multisine experiments
% System:   Conventional motor-bench with flexible coupling
% Author:   Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('private/MultisineTypeB.mat');     % random odd-odd multisine experiment

% reproduce input
% harm.fs = fs;
% harm.df = df;
% harm.fl = fl;
% harm.fh = fh;
% harm.fr = fr;
% Hampl = repmat(tf(1),[1,1]); % flat spectrum
% options.itp = itp;
% options.ctp = ctp;
% options.dtp = stp;
% options.gtp = 'l';      % grid type: l=linear/q=quasi-logarithmic
% ms = multisine(harm, Hampl, options);


% Time treatment: remove transients/offsets/trends
trans = 1;                      % number of transient periods
trend = 0;                      % period trend removal flag
input = [iq_adx];               % input data to motor bench
output = [theta_mx,-theta_my];  % output data of motor bench
[x,time] = pretreat(input,nrofs,fs,trans,trend);
[y,time] = pretreat(output,nrofs,fs,trans,trend);

%% STEP 3: NON-LINEARITY DETECTION
% quantify non-linear distortion level in data
[Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,fs,fl,fh,df);

figure
subplot(211), semilogx(freql,dbm(Yl(:,1)),freqe,dbm(Ye(:,1)),...
                       freqo,dbm(Yo(:,1)),freqn,dbm(Yn(:,1)));
    ylabel('magnitude [dB]'), xlim([fl,fh])
    legend('Ylin','Yeven','Yodd','Ynoise')
subplot(212), semilogx(freql,dbm(Yl(:,2)),freqe,dbm(Ye(:,2)),...
                       freqo,dbm(Yo(:,2)),freqn,dbm(Yn(:,2)));
    ylabel('magnitude [dB]'), xlim([fl,fh])
    
% NOTE: important low-frequent odd-distortions present in data
%       probably due to non-linear rolling friction in the bearings 