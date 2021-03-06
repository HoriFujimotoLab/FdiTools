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
    xlabel('freq [Hz]'); ylabel('amplitude [dB]');

%% EXPERIMENT
input = ss.x*1.6726;
load('private/20160829_ident'); % load benchmark model
inputnoize = 0.01; % amp of input noise 
input = input + inputnoize*randn(size(input));
Ts = 1/harm.fs;
t = 0:Ts:Ts*(length(input)-1);
output = lsim(mdl.Pv(1,1),input,t);
outputnoize = 0.001; % amp of output noise 
output = output + outputnoize*randn(size(output));

%% STEP 2: NonparametricFRF
% remove transient periods, offsets and trends
input = detrend(input,0); output = detrend(output,0);
[txy,freq] = tfestimate(input,output,rectwin(harm.df*harm.fs),0,harm.df*harm.fs,harm.fs);
[cxy,freq] = mscohere(input,output,rectwin(harm.df*harm.fs),0,harm.df*harm.fs,harm.fs);
Pfrd = frd(txy,freq,'FrequencyUnit','Hz');
figure; semilogx(freq,cxy); title('coherence');

% require system identification toolbox
Pfrd2 = fselect(Pfrd,harm.df,harm.fh);
opt = tfestOptions('WeightingFilter',cxy.*freq);
Pest = tfest(Pfrd2,7,4);
% data = iddata(output,input,1/fs); 
% Pest = tfest(data,7,4); % require system identification toolbox
figure; bode(Pfrd,Pest,mdl.Pv(1,1)); xlim([1,1000]); % require system identification toolbox
legend('estimated FRF','fitted by tfest','TRUE');