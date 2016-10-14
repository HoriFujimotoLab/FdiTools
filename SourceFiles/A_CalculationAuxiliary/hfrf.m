function [FRF] = hfrf(Hm,freq)
%HFRF - get frequency response matrix of model.
%   [FRF] = hfrf(Hm,freq)
% Hm    : structure with models to evaluate.
% FRF   : structure containing calculated frf.
% author: Thomas Beauduin, KULeuven, PMA, 2014
%%%%%
nrofo = size(Hm,1);         % number of output
nrofi = size(Hm,2);         % number of input
nrofh = nrofi*nrofo;        % number of transfer functions
nroff = length(freq(:));    % number of freq lines

FRF = zeros(nroff,nrofh);
for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    FRF(:,h) = squeeze(freqresp(Hm(o,i),freq*2*pi));
end

end

