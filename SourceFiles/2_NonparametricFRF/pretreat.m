function [y,time] = pretreat(x,nrofs,varargin)
%PRETREAT - pretreatment of data vectors. (MIMO)
% [y,time] = pretreat(x,nrofs,fs)
% x     : untreated periodic measurement vector
% nrofs : number of samples per signal period
% fs    : sampling frequency
% nroft : number of transient periods
% trend : flag for trend removal {0,1}
% y     : transient, offset and trend removed data
% Author: Thomas Beauduin, KULeuven, 2014
%%%%%
switch nargin-2                         % Input processing
    case 0, fs = 1; nroft = 1; trend = 0;
    case 1, fs = varargin{1}; nroft = 1; trend = 0;
    case 2, fs = varargin{1}; nroft = varargin{2}; trend = 0;
    case 3, fs = varargin{1}; nroft = varargin{2}; trend = varargin{3};
end

nrofr = size(x,1);                      % number of data points (rows)
nrofc = size(x,2);                      % number of data vectors (columns)
nrofp = ceil(nrofr/nrofs)-nroft;        % number of periods
y = x(1:nrofs*nrofp,:);

for k=1:nrofc
    u = x(nroft*nrofs+1:end,k);         % Transient Removal
    D = reshape(u,nrofs,nrofp);         % Offset Removal
    Da = mean(D,1);
    for m=1:nrofp
        for n=1:nrofs
            D(n,m) = D(n,m) - Da(m);
        end
    end
    if trend > 0                        % Trend Removal
        for m=1:nrofp
            D(:,m) = detrend(D(:,m));
        end
    end
    y(:,k) = D(:);
end

time = (0:(length(y(:,1))-1))'/fs;
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