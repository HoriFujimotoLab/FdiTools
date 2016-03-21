function [cf,cferr]=crestfactor(X,N,Fe)
% CRESTFACTOR - Compute the signal Crest factors.
%
% X         : data fft matrix with signals column by column
% N         : number of time samples to be used
% Fe        : set of the effective harmonic numbers
% cf        : lower and upper crest factor values
% cferr     : lower and upper crest factor errors
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
% see also TIMEFACTOR, SCHROEDER

if nargin<3
   if nargin<2, x=f2t(X); else x=f2t(X,N); end
   cf = lpnorm(x,inf)./effval(X);
else
   if isempty(N), x=f2t(X); else x=f2t(X,N); end
   cf = lpnorm(x,inf)./effval(X,Fe);
end
if nargin<2, [N,nCol]=size(x); end
[nRow,nCol]=size(X);
if nCol==1, dummy=X; else dummy=max(abs(X)')'; end
Ft = cumsum(ones(size(dummy)));
Ft = Ft(dummy>max(dummy)*1e-6);
clear dummy
if nargin<3, Fe=Ft; end
kmax=max(max(Ft))-1;
if N>pi*kmax
   cferr=cf/(1-pi*kmax/N);
end

