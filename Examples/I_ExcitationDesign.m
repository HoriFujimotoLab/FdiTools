%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCITATION DESIGN:
% ------------------
% Descr.:   example of several excitation signals design
%           to demonstrate the design parameters/options
% Author:   Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% EXAMPLE 1: MULTISINE
% Harmonics Parameters:
harm.fs = 1000;         % sampling frequency
harm.df = 1;            % frequency resolution
harm.fl = 1;            % lowest frequency
harm.fh = 100;          % highest frequency
harm.fr = 1.02;         % frequency log ratio
% Design Options:
options.itp = 'r';      % init phase type:  s=schroeder/r=random
options.ctp = 'c';      % compression type: c=comp/n=no_comp
options.dtp = 'o';      % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
options.gtp = 'l';      % grid type: l=linear/q=quasi-logarithmic
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
pause

% NOTES:
% CTP: Compression type algorithm optimizes phase for 'min crest factor' 
%      to improve S/N of measurement.
% ITP: Random initial phase creates different signals in time domain
%      with identical frequency domain.
% DTP: Odd signal type used for non-linear distortion analysis.
% GTP: Quasi-logarithmic freq grid used for overview measurements

%% EXAMPLE 2: PRBS
% Experiment Parameters:
fs = 1000;              % sampling frequency   [Hz]
df = 1;                 % frequency resolution [Hz]
fl = 10;                % min excitation freq. [Hz]
fh = 100;               % max excitation freq. [Hz]
nrofs = fs/df; time = (0:1/fs:1/df-1/fs);

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
% Harmonic Parameters
harm.fs = 2000;             % sampling frequency   [Hz]
harm.df = 0.05;             % frequency resolution [Hz]
harm.fl = 0.1;              % lowest frequency     [Hz]
harm.fh = 100;              % highest frequency    [Hz]

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
