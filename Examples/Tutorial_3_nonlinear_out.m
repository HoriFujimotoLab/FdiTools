%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TUTORIAL OF FDI by QLOG:
% ----------------------
% Descr.:   Tutorial of FdiTools by qlog excitation.
% System:   High-precision positioning stages with two encoders.
% Author:   Wataru Ohnishi, The University of Tokyo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% STEP 1: ExcitationDesign
harm.fs = 10000;         % sampling frequency
harm.df = 1;            % frequency resolution
harm.fl = 1;            % lowest frequency
harm.fh = 1000;          % highest frequency
harm.fr = 1.02;         % frequency log ratio
% Design Options:
options.itp = 'r';      % init phase type:  s=schroeder/r=random
options.ctp = 'c';      % compression type: c=comp/n=no_comp
options.dtp = 'O';      % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
options.gtp = 'f';      % grid type: l=linear/q=quasi-logarithmic
% Ampliude spectrum:
nrofi = 1;              % Define number of inputs
Hampl = repmat(tf(1),[1,nrofi]); % flat spectrum
ms = multisine(harm, Hampl, options);

hfig=figure; sub = 0;
for ii = 1:nrofi
    for jj = 1:nrofi
        sub = sub+1;
        subplot(nrofi, nrofi, sub);
        plot(ms.time,squeeze(ms.x(ii,jj,:)))
        title(strcat('CF =',num2str(ms.cf(ii,jj))))
        old = axis; axis([0, ms.time(end), old(3:4)])
    end
end
hfig=figure; sub=0;
for ii = 1:nrofi
    for jj = 1:nrofi
        sub = sub+1;
        subplot(nrofi, nrofi, sub);
        semilogx(ms.freq,dbm(squeeze(ms.X(ii,jj,:))),'+')
        old = axis; axis([0, length(ms.x), old(3:4)])
    end
end

%% EXPERIMENT
load('private/20160829_ident'); % load benchmark model
G = mdl.Pv(1,1);
rng default
nrofs = length(ms.x(1,1,:));
input = squeeze(ms.x(1,1,:));
nrofp = 20; % number of period of periodic excitation
input = repmat(input,[nrofp,1]);
input_noise = 0*randn(size(input)); % amp of input noise 
Ts = 1/ms.harm.fs;
t = 0:Ts:Ts*(length(input)-1);
odd_nl_amp = 0.1;
even_nl_amp = 1;
output_noise = 0.0001*randn(size(input)); % amp of output noise 

% exp1
input_amp1 = 1; % amp of input 
input1 = input_amp1*input; 
input1_noise = input1 + input_noise;
in.time = t; in.signals.values = input1_noise;
temp = sim('model_nl_out');
output1_noise = temp.out + output_noise;

% exp2
input_amp2 = 0.1; % amp of input 
input2 = input_amp2*input; 
input2_noise = input2 + input_noise;
in.time = t; in.signals.values = input2_noise;
temp = sim('model_nl_out');
output2_noise = temp.out + output_noise;

% exp1
input_amp3 = 10; % amp of input 
input3 = input_amp3*input; 
input3_noise = input3 + input_noise;
in.time = t; in.signals.values = input3_noise;
temp = sim('model_nl_out');
output3_noise = temp.out + output_noise;

% noise-free condition
output_noise_free = lsim(mdl.Pv(1,1),input1,t);
%% STEP 3: NON-LINEARITY DETECTION
% input 1
% remove transient periods, offsets and trends
trans = 2;                      % number of transient periods
trend = 0;                      % period trend removal flag
[x,time] = pretreat(input1,nrofs,harm.fs,trans,trend);
[y,time] = pretreat(output1_noise,nrofs,harm.fs,trans,trend);

% quantify non-linear distortion level in data
[Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,harm.fs,harm.fl,harm.fh,harm.df);

figure;
semilogx(freql,dbm(Yl(:,1)),'*',freqe,dbm(Ye(:,1)),'*',...
                       freqo,dbm(Yo(:,1)),'*',freqn,dbm(Yn(:,1)),'*');
ylabel('Magnitude [dB]'), xlim([harm.fl,harm.fh])
legend('Ylin','Yeven','Yodd','Ynoise')
title('with noise and nl with input amp 1');

% input 0.1
% remove transient periods, offsets and trends
trans = 2;                      % number of transient periods
trend = 0;                      % period trend removal flag
[x,time] = pretreat(input2,nrofs,harm.fs,trans,trend);
[y,time] = pretreat(output2_noise,nrofs,harm.fs,trans,trend);

% quantify non-linear distortion level in data
[Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,harm.fs,harm.fl,harm.fh,harm.df);

figure;
semilogx(freql,dbm(Yl(:,1)),'*',freqe,dbm(Ye(:,1)),'*',...
                       freqo,dbm(Yo(:,1)),'*',freqn,dbm(Yn(:,1)),'*');
ylabel('Magnitude [dB]'), xlim([harm.fl,harm.fh])
legend('Ylin','Yeven','Yodd','Ynoise')
title('with noise and nl with input amp 0.1');

% input 10
% remove transient periods, offsets and trends
trans = 2;                      % number of transient periods
trend = 0;                      % period trend removal flag
[x,time] = pretreat(input3,nrofs,harm.fs,trans,trend);
[y,time] = pretreat(output3_noise,nrofs,harm.fs,trans,trend);

% quantify non-linear distortion level in data
[Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,harm.fs,harm.fl,harm.fh,harm.df);

figure;
semilogx(freql,dbm(Yl(:,1)),'*',freqe,dbm(Ye(:,1)),'*',...
                       freqo,dbm(Yo(:,1)),'*',freqn,dbm(Yn(:,1)),'*');
ylabel('Magnitude [dB]'), xlim([harm.fl,harm.fh])
legend('Ylin','Yeven','Yodd','Ynoise')
title('with noise and nl with input amp 10');

% noise-free condition
trans = 2;                      % number of transient periods
trend = 0;                      % period trend removal flag
[x,time] = pretreat(input,nrofs,harm.fs,trans,trend);
[y,time] = pretreat(output_noise_free,nrofs,harm.fs,trans,trend);

[Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,harm.fs,harm.fl,harm.fh,harm.df);

figure;
semilogx(freql,dbm(Yl(:,1)),'*',freqe,dbm(Ye(:,1)),'*',...
                       freqo,dbm(Yo(:,1)),'*',freqn,dbm(Yn(:,1)),'*');
ylabel('Magnitude [dB]'), xlim([harm.fl,harm.fh])
legend('Ylin','Yeven','Yodd','Ynoise');
title('without noise and nl');

