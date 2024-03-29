%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TUTORIAL OF FDI by QLOG:
% ----------------------
% Descr.:   Tutorial of FdiTools by qlog excitation.
% System:   High-precision positioning stages with two encoders.
% Author:   Wataru Ohnishi, The University of Tokyo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
rng('default');

% Experiment settings
load('private/20160829_ident'); % load benchmark model
P0 = mdl.Pv(1,1);
inputnoize = 0.01; % amp of input noise
outputnoize = 0.003; % amp of output noise

%% First experiment with wideband excitation (quasi log)
% STEP 1: ExcitationDesign
harm.fs = 10000;         % sampling frequency
harm.df = 1;            % frequency resolution
harm.fl = 1;            % lowest frequency
harm.fh = 1000;          % highest frequency
harm.fr = 1.02;         % frequency log ratio
% Design Options:
options.itp = 'r';      % init phase type:  s=schroeder/r=random
options.ctp = 'c';      % compression type: c=comp/n=no_comp
options.dtp = 'f';      % signal type:      f=full/ O=odd-odd
%                   o=odd / O2=special odd-odd
options.gtp = 'q';      % grid type: l=linear/q=quasi-logarithmic
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

% EXPERIMENT
nrofs = length(ms.x(1,1,:));
input = squeeze(ms.x(1,1,:))/(abs(max(squeeze(ms.x(1,1,:)))));
nrofp = 5; % number of period of periodic excitation
input = repmat(input,[nrofp,1]);
input = input + inputnoize*randn(size(input));
Ts = 1/ms.harm.fs; t = 0:Ts:Ts*(length(input)-1);
output = lsim(P0,input,t); % EXPERIMENT
output = output + outputnoize*randn(size(output));

% STEP 2: NonparametricFRF
% remove transient periods, offsets and trends
trans = 2;                      % number of transient periods
trend = 0;                      % period trend removal flag
r0 = (1:ms.nrofs);                 % time visualization range
rn = (ms.nrofs*2+1:(2+1)*ms.nrofs);   % data visualization range
[x,time] = pretreat(input,nrofs,harm.fs,trans,trend);
[y,time] = pretreat(output,nrofs,harm.fs,trans,trend);

Pest = time2frf_ml(x,y,ms);
bode_fdi({P0,Pest});
legend('true','estimated frd','sGhat');

disp('----------------------------------------')
disp('First experiment with wideband excitation')
disp('----------------------------------------')

Pest_qlog = Pest;
Pest_qlog_del = fdel_fdi(Pest_qlog,100,1000);

pause


%% Second experiment via inverse S/N ratio
% STEP 1: ExcitationDesign
harm.fs = 10000;         % sampling frequency
harm.df = 1;            % frequency resolution
harm.fl = 101;            % lowest frequency
harm.fh = 296;;          % highest frequency
harm.fr = 1.02;         % frequency log ratio
% Design Options:
options.itp = 'r';      % init phase type:  s=schroeder/r=random
options.ctp = 'c';      % compression type: c=comp/n=no_comp
options.dtp = 'f';      % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
options.gtp = 'q';      % grid type: l=linear/q=quasi-logarithmic
% Ampliude spectrum:
nrofi = 1;              % Define number of inputs
sGhat_inv = frd((abs(squeeze(Pest_qlog.UserData.sGhat))./abs(squeeze(Pest_qlog.resp))),Pest_qlog.freq,'FrequencyUnit','Hz');
ms = multisine(harm, sGhat_inv, options);

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

% EXPERIMENT
nrofs = length(ms.x(1,1,:));
input = squeeze(ms.x(1,1,:));
nrofp = 20; % number of period of periodic excitation
input = repmat(input,[nrofp,1]);
input = input + inputnoize*randn(size(input));
Ts = 1/ms.harm.fs; t = 0:Ts:Ts*(length(input)-1);
output = lsim(P0,input,t); % EXPERIMENT
output = output + outputnoize*randn(size(output));

% STEP 2: NonparametricFRF
% remove transient periods, offsets and trends
trans = 2;                      % number of transient periods
trend = 0;                      % period trend removal flag
r0 = (1:ms.nrofs);                 % time visualization range
rn = (ms.nrofs*2+1:(2+1)*ms.nrofs);   % data visualization range
[x,time] = pretreat(input,nrofs,harm.fs,trans,trend);
[y,time] = pretreat(output,nrofs,harm.fs,trans,trend);

figure
subplot(211), plot(time(r0),x(rn,:))
    title('input data'), legend('input force'), ylim([-2,2])
subplot(212), plot(time(r0),y(rn,:))
    title('output data'), legend('stage position')
    xlabel('time [s]')

Pest = time2frf_ml(x,y,ms);
bode_fdi({P0,Pest(1,1)},[Pest.freq,Pest.UserData.FRFn(:,1)]);
legend('true','estimated frd','noise');
Pest_lin_100_300 = Pest;

disp('----------------------------------------')
disp('Second experiment via inverse S/N ratio')
disp('----------------------------------------')

pause

%% Third experiment around high-frequency range
% STEP 1: ExcitationDesign
harm.fs = 10000;         % sampling frequency
harm.df = 1;            % frequency resolution
harm.fl = 302;            % lowest frequency
harm.fh = 1e3;          % highest frequency
harm.fr = 1.02;%1.005;         % frequency log ratio
% Design Options:
options.itp = 'r';      % init phase type:  s=schroeder/r=random
options.ctp = 'c';      % compression type: c=comp/n=no_comp
options.dtp = 'f';      % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
options.gtp = 'q';      % grid type: l=linear/q=quasi-logarithmic
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

% EXPERIMENT
nrofs = length(ms.x(1,1,:));
input = squeeze(ms.x(1,1,:));
nrofp = 20; % number of period of periodic excitation
input = repmat(input,[nrofp,1]);
input = input + inputnoize*randn(size(input));
Ts = 1/ms.harm.fs; t = 0:Ts:Ts*(length(input)-1);
output = lsim(P0,input,t); % EXPERIMENT
output = output + outputnoize*randn(size(output));

% STEP 2: NonparametricFRF
% remove transient periods, offsets and trends
trans = 2;                      % number of transient periods
trend = 0;                      % period trend removal flag
r0 = (1:ms.nrofs);                 % time visualization range
rn = (ms.nrofs*2+1:(2+1)*ms.nrofs);   % data visualization range
[x,time] = pretreat(input,nrofs,harm.fs,trans,trend);
[y,time] = pretreat(output,nrofs,harm.fs,trans,trend);

figure
subplot(211), plot(time(r0),x(rn,:))
    title('input data'), legend('input force'), ylim([-2,2])
subplot(212), plot(time(r0),y(rn,:))
    title('output data'), legend('stage position')
    xlabel('time [s]')

Pest = time2frf_ml(x,y,ms);
bode_fdi({P0,Pest(1,1)},[Pest.freq,Pest.UserData.FRFn(:,1)]);
legend('true','estimated frd','noise');
Pest_qlog_301_1000 = Pest;

disp('----------------------------------------')
disp('Third experiment around high-frequency range')
disp('----------------------------------------')

pause
%% Connect all experiment
bode_fdi({P0,Pest_qlog(1,1)},[Pest_qlog.freq,Pest_qlog.UserData.sGhat(:,1)]);
legend('true','estimated frd','noise');
title('Single experiment');

Pest = fcat_fdi(Pest_qlog_del,Pest_lin_100_300,Pest_qlog_301_1000);
bode_fdi({P0,Pest(1,1)},[Pest.freq,Pest.UserData.sGhat(:,1)]);
legend('true','estimated frd','noise');
title('Iterative experiment');

EstErr_original = 1-P0/Pest_qlog;
EstErr_after = 1-P0/Pest;
figure; bodemag(EstErr_original,'r',EstErr_after,'b');
title('Estimation error');
legend('Single experiment','Iterative experiment','Location','northwest');

disp('----------------------------------------')
disp('Connect three expeiment via sGhat')
disp('----------------------------------------')

pause
%% STEP 4: Parametric estimation
% deterministric/stochastic estimation with non-parametric noise model
n=7;                        % model order of denominator polynomial
mh=[4;]; ml=[0;];     % model orders of numerator polynomial
relvar=1e-10;               % relative variation of costfunction (stop)
iter=5e2;                   % maximum number of iterations (stop)
GN = 0;                     % Levenberg-Marquardt optimization
cORd = 'c';                 % continuous model identifaction
FRF_W = ones(size(squeeze(Pest.resp)));   % least squares weighting function
relax = 1;                  % relaxation factor for btls estimation
freq = Pest.freq;

% Deterministic methods: non-linear least squares
[SYS.nls,SYS.wls] = nlsfdi(Pest,FRF_W,n,mh,ml,iter,relvar,GN,cORd);
    FRF.nls = hfrf(SYS.nls,freq);
    FRF.wls = hfrf(SYS.wls,freq);
% Stochastic methods: maximum likelihood estimation
[SYS.ml,SYS.ls]  = mlfdi(Pest,n,mh,ml,iter,relvar,GN,cORd);
    FRF.ml = hfrf(SYS.ml,freq);
    FRF.ls = hfrf(SYS.ls,freq);

% Stochastic methods: bootstrapped total least squares
[SYS.btls,SYS.gtls] = btlsfdi(Pest,n,mh,ml,relax,iter,relvar,cORd);
    FRF.btls = hfrf(SYS.btls,freq);
    FRF.gtls = hfrf(SYS.gtls,freq);

bode_fdi({P0,Pest(1,1),SYS.wls,SYS.nls,SYS.ls},[Pest.freq,Pest.UserData.FRFn(:,1)]);
legend('TRUE','FRF','WLS','NLS','LS','FRFn');

bode_fdi({P0,Pest(1,1),SYS.ml,SYS.btls,SYS.gtls},[Pest.freq,Pest.UserData.sGhat(:,1)]);
legend('TRUE','FRF','MLE','BTLS','GTLS','sGhat')

% best estimator
bode_fdi({P0,Pest(1,1),SYS.btls},[Pest.freq,Pest.UserData.sGhat(:,1)]);
legend('TRUE','FRF','BTLS','sGhat')

disp('----------------------------------------')
disp('Parametric estimation')
disp('----------------------------------------')