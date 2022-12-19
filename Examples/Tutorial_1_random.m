%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TUTORIAL OF RANDOM NOISE EXCITATION
% ----------------------
% Descr.:   Tutorial of random noise excitation. (conventional)
% System:   High-precision positioning stages with two encoders.
% Author:   Wataru Ohnishi, The University of Tokyo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% STEP 1: ExcitationDesign
fs = 10000;         % sampling frequency
texp = 5;           % experiment time
amp = 1.6726;       % max amplitude of multisine for comparison
rng('default');     % for reproducibility
input = rand(fs*texp,1)*2*amp -amp;

%% EXPERIMENT
load('private/20160829_ident'); % load benchmark model
P0 = mdl.Pv(1,1);
inputnoize = 0.01; % amp of input noise 
input = input + inputnoize*randn(size(input));
Ts = 1/fs;
t = 0:Ts:Ts*(length(input)-1);
input_noise_amp = 0.01; % amp of input noise 
input_noize = input + input_noise_amp*randn(size(input));
output = lsim(P0,input_noize,t);
output_noise_amp = 0.001; % amp of output noise 
output_noize = output + output_noise_amp*randn(size(output));


%% STEP 2: NonparametricFRF
% remove transient periods, offsets and trends
input = detrend(input,0);output_noize = detrend(output_noize,0);
[txy,freq] = tfestimate(input,output_noize,[],[],[],fs);
[cxy,freq] = mscohere(input,output_noize,[],[],[],fs);
Pfrd = frd(txy,freq,'FrequencyUnit','Hz');
figure; semilogx(freq,cxy); title('coherence');

% require system identification toolbox
if exist('tfestOptions')
    opt = tfestOptions('WeightingFilter',cxy.*freq);
    Pest = tfest(Pfrd,7,4);
    % data = iddata(output,input,1/fs); 
    % Pest = tfest(data,7,4); % require system identification toolbox
    bop = bodeoptions('cstprefs');
    bop.PhaseWrapping = 'on';
    figure; bode(Pfrd,Pest,P0,bop); xlim([1,1000]); % require system identification toolbox
    legend('estimated FRF','fitted by tfest','TRUE');
end