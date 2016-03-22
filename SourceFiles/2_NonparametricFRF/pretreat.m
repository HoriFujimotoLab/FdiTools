function [y,time] = pretreat(x,nrofs,varargin)
%PRETREAT - pretreatment of data vectors for estimation.
% [y,time] = pretreat(x,nrofs,fs)
% x     : untreated periodic measurement vector
% nrofs : number of samples per signal period
% fs    : sampling frequency
% nroft : number of transient periods
% trend : flag for trend removal {0,1}
% y     : transient, offset and trend removed data
% Author: Thomas Beauduin, KULeuven, 2014
%%%%%

% 0. Input processing:
switch nargin-2
    case 0, fs = 1; nroft = 1; trend = 0;
    case 1, fs = varargin{1}; nroft = 1; trend = 0;
    case 2, fs = varargin{1}; nroft = varargin{2}; trend = 0;
    case 3, fs = varargin{1}; nroft = varargin{2}; trend = varargin{3};
end

% 1. Transient Removal:
x = x(nroft*nrofs+1:end);
nrofp = ceil(length(x)/nrofs);

% 2. Offset Removal:
P = reshape(x,nrofs,nrofp);
Pa = mean(P,1);
for k=1:nrofp
    for l=1:nrofs
        P(l,k)=P(l,k)-Pa(k);
    end
end
y = P(:);
time = (0:(length(y)-1))'/fs;

% 3. Trend Removal:
if trend > 0
    P = reshape(y,nrofs,nrofp);
    for k=1:nrofp
        P(:,k)=detrend(P(:,k));
    end
    y = P(:);
end
end

% 4. References:
% McCormack, A. S., J. O. Flower, and K. R. Godfrey (1994). 
% The suppression of drift and transient effects for 
% frequecy-domain identification. 
% IEEE Trans. Instrum. and Meas., vol. 43, no. 2, pp.232-237.

% Peirlinckx, L., P. Guillaume, and R. Pintelon (1996). 
% Accurate and fast estimation of the Fourier coefficients of 
% periodic signals disturbed by a trend. 
% IEEE Trans. Instrum. and Meas., vol. 45, no.1, pp. 5-11.