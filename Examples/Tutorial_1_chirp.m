%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TUTORIAL OF RANDOM NOISE EXCITATION
% ----------------------
% Descr.:   Tutorial of random noise excitation. (conventional)
% System:   High-precision positioning stages with two encoders.
% Author:   Wataru Ohnishi, The University of Tokyo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% STEP 1: ExcitationDesign
harm.fs = 10000;            % sampling frequency   [Hz]
harm.df = 1;             % frequency resolution [Hz]
harm.fl = 1;              % lowest frequency     [Hz]
harm.fh = 1000;              % highest frequency    [Hz]

% Design Options:
options.type = 'lin';       % Sweep type: lin=linear/qdr=quadratic
%             log=logarithmic
ss = sweptsine(harm,options);

figure
subplot(211); plot(ss.time,ss.x);
title('swept-sine: time domain');
xlabel('time [s]'); ylabel('amplitude [-]');
subplot(212); semilogx(ss.freq,dbm(ss.X));
title('Swept-sine: freq domain');
xlabel('Frequency [Hz]'); ylabel('amplitude [dB]');

%% EXPERIMENT
nrofs = length(ss.x); % number of sample per period
input = ss.x*1.6726;
nrofp = 10; % number of period of periodic excitation
input = repmat(input,[nrofp,1]);
load('private/20160829_ident'); % load benchmark model
input_noise_amp = 0.01; % amp of input noise
input_noize = input + input_noise_amp*randn(size(input));
Ts = 1/harm.fs;
t = 0:Ts:Ts*(length(input)-1);
output = lsim(mdl.Pv(1,1),input_noize,t);
output_noise_amp = 0.001; % amp of output noise
output_noize = output + output_noise_amp*randn(size(output));

%% STEP 2: NonparametricFRF
% remove transient periods, offsets and trends
trans = 1;                      % number of transient periods
trend = 0;                      % period trend removal flag
[input_pretreat,time] = pretreat(input,nrofs,harm.fs,trans,trend);
[output_pretreat,time] = pretreat(output_noize,nrofs,harm.fs,trans,trend);

[txy,freq] = tfestimate(input_pretreat,output_pretreat,rectwin(harm.df*harm.fs),0,harm.df*harm.fs,harm.fs);
[cxy,freq] = mscohere(input_pretreat,output_pretreat,rectwin(harm.df*harm.fs),0,harm.df*harm.fs,harm.fs);
Pfrd = frd(txy,freq,'FrequencyUnit','Hz');
figure; semilogx(freq,cxy);
xlabel('Frequency [Hz]'); ylabel('Coherence [-]');

% require system identification toolbox
if exist('tfestOptions')
    Pfrd2 = fselect(Pfrd,harm.df,harm.fh);
    opt = tfestOptions('WeightingFilter',cxy.*freq);
    Pest = tfest(Pfrd2,7,4);
    bop = bodeoptions('cstprefs');
    bop.PhaseWrapping = 'on';
    figure; bode(Pfrd,Pest,mdl.Pv(1,1),bop); xlim([1,1000]);
    legend('estimated FRF','fitted by tfest','TRUE');
end