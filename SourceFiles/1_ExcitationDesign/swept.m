function [x,time,X,freq] = swept(fs,fl,fh,df)
%SWEPT - Swept-Sine excitation signal generation.
%   [x,time,X,freq] = swept(fs,fl,fh,df)
% fs, df    : sampling frequency & resolution [Hz]
% fl,fh     : lower/upper frequency range limit [Hz]
% x, time   : swept sine signal time domain [s]
% Algorithm : x(t) = Asin((at+b)t)   for  0<=t<T
% Author    : Thomas Beauduin, KULeuven, PMA, 2014
%%%%%

nrofs = fs/df;
time = (0:(nrofs-1))'/fs;
k1 = ceil(fl/df); k2 = floor(fh/df);

a = pi*(k2-k1)*df^2;
b = 2*pi*k1*df;
x = sin((a*time+b).*time);

X=t2f(x,nrofs);
freq=fs*(0:1:nrofs/2-1)'/nrofs; 