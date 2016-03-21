function [x,time,X,freq] = swept(fs,fl,fh,df)
%SWEPT - Swept-Sine excitation signal generation.
%
% fs        : sample frequency [Hz]
% fl        : lower limit of selected frequency range [Hz]
% fh        : upper limit of selected frequency range [Hz]
% df        : distance between spectral lines
% x         : swept sine signal
% time      : time vector [s]
% Algorithm : x(t) = Asin((at+b)t)   for  0<=t<T
% Author    : Thomas Beauduin, KULeuven, 2014
%%%%%

nrofs = fs/df;
time = [0:(nrofs-1)]'/fs;
k1 = ceil(fl/df); k2 = floor(fh/df);

a = pi*(k2-k1)*df^2;
b = 2*pi*k1*df;
x = sin((a*time+b).*time);

X=t2f(x,nrofs);
freq=fs*(0:1:nrofs/2-1)'/nrofs; 