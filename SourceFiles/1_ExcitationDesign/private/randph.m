function X = randph(X) 
%RANDPH - random multisine phase design [-pi,pi]
%   X = randph(X)
% X (in)   : ampl. spectrum of signals column by column
% X (out)  : spectrum with random phases 
% Author   : Thomas Beauduin, KULeuven, 2014
%%%%%
ampl = abs(X);
rand('seed',sum(100*clock)); 
phase = (2*rand(size(ampl))-1)*pi; 
X = ampl.*exp(1i*phase); 
 
end