%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCITATION DESIGN:
% ------------------
% Descr.:   example of several excitation signals design
%           to demonstrate the design parameters/options
% Author:   Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% Experiment Parameters:
fs = 1000;              % sampling frequency   [Hz]
df = 1;                 % frequency resolution [Hz]
fl = 10;                % min excitation freq. [Hz]
fh = 100;               % max excitation freq. [Hz]
nrofs = fs/df; time = (0:1/fs:1/df-1/fs);

%% EXAMPLE 1: MULTISINE
% Type A: full compressed schroeder-phase multisine:
itp = 's';              % init phase type:  s=schroeder/r=random
ctp = 'c';              % compression type: c=comp/n=no_comp
stp = 'f';              % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
Bn=1; An=1;             % amplitude spectrum    [-]
[x1,Xs1,freqs1,Xt1,freqt1] = msin(fs,df,fl,fh,itp,ctp,stp,Bn,An);

% Type B: odd-odd uncompressed random-phase multisine:
itp = 'r';              % init phase type:  s=schroeder/r=random
ctp = 'n';              % compression type: c=comp/n=no_comp
stp = 'O';              % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
Bn=1; An=1;             % amplitude spectrum (opt) [-]
[x2,Xs2,freqs2,Xt2,freqt2] = msin(fs,df,fl,fh,itp,ctp,stp,Bn,An);

figure
subplot(211); plot(time,[x1(1:nrofs),x2(1:nrofs)]);
    title('multisine: time domain'); 
    xlabel('time [s]'); ylabel('amplitude [-]');
legend('Type A','Type B')
subplot(212); semilogx(freqt1,[dbm(Xt1),dbm(Xt2)]);
    title('multisine: freq domain'); 
    xlabel('freq [Hz]'); ylabel('amplitude [dB]');
pause

% NOTES:
% CTP: Compression type algorithm optimizes phase for 'min crest factor' 
%      to improve S/N of measurement.
% ITP: Random initial phase creates different signals in time domain
%      with identical frequency domain.
% STP: Odd signal type used for non-linear distortion analysis.
%
%% EXAMPLE 2: PRBS 
% Pseudo-Random-Binary-Sequency excitation signal
log2N = 21;             % shift reg length     [-]
bitno = fs/df;          % number of bits       [-]
[x,X,freq,nextstnum] = prbs(fs,log2N,bitno);

figure
subplot(211); plot(time,x);
    title('prbs: time domain'),xlabel('time [s]')
    ylabel('amplitude [-]'),ylim([-2,2])
subplot(212); semilogx(freq,dbm(X));
    title('freq domain signal'),xlabel('frequency [Hz]')
    ylabel('amplitude [dB]'),xlim([1,freq(end)])
pause

%% EXAMPLE 3: SWEPT-SINE
% Swept-Sine excitation signal generation
[x,time,X,freq] = swept(fs,fl,fh,df);

figure
subplot(211); plot(time,x(1:nrofs));
    title('swept-sine: time domain'); 
    xlabel('time [s]'); ylabel('amplitude [-]');
subplot(212); semilogx(freq,dbm(X));
    title('Swept-sine: freq domain'); 
    xlabel('freq [Hz]'); ylabel('amplitude [dB]');
