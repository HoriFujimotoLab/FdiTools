function [output] = sweptsine(h, options)
%SWEPTSINE - Swept-Sine excitation signal generation.
%
% HARMONICS = parameter set containing excitation harmonics
%   <>.fs     : Sampling frequency in Hz
%   <>.fl/.fh : Lowest/Highest frequency lines
%   <>.df     : Spectral lines density
% OPTIONS = parameter set containing design options
%   <>.type   : sweep - 'lin' linear/'qdr' quadratic/'log' logarithmic
% OUTPUT = time/frequency information of designed signal
%   <>.x/time : time domain data of multisine
%   <>.X/freq : frequency domain data of multisine
% Author      : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%

% Calculation of signal parameters
nrofs = ceil(h.fs/h.df);                % Number of samples/period
time = (0:(nrofs-1))'/h.fs;             % Time vector of data
freq = h.fs*(0:1:nrofs/2-1)'/nrofs;     % Frequency vector

% Calculation of sweptsine signal type
switch options.type
    case {'lin','linear'},
        x = sin( pi*(h.fh-h.fl)*h.df.*(time.^(2))+2*pi*h.fl.*time );
    case {'qdr','quadratic'},
        x = sin( 2/3*pi*(h.fh-h.fl)*h.df^2.*(time.^(3))+2*pi*h.fl.*time );
    case {'log','logarithmic'},
        x = sin( 2*pi/h.df/log(h.fh/h.fl)* (h.fl*(h.fh/h.fl).^(time*h.df)-h.fl) );
end

% Calculation of frequency domain
X = fft(x,nrofs);
X = [X(1,:);2*X(2:floor(nrofs/2),:)]/nrofs;

% Ouput structure creation
output.x = x; output.time = time;
output.X = X; output.freq = freq; 

end


% Consider adding crest-factor and time-factor calculation