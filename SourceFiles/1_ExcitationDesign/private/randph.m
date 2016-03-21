function X=randph(X) 
%RANDPH - random multisine phase design [-pi,pi]
%
% X (in)   : abs(X) ampl. spectrum of the signals column by column. 
% X (out)  : spectrum with random phases 
% Algorithm: phase.n = phase.1 - 2p.[S{k=1,n-1} (n-k)Ak^2]
% Author   : Thomas Beauduin, KULeuven, 2014
% 
ampl=abs(X);
rand('seed',sum(100*clock)); 
phase = (2*rand(size(ampl))-1)*pi; 
X=ampl.*exp(1i*phase); 
 
