function eff = effval(X,Fe)
% EFFVAL - Calculates effective value of signal.
%
% X         : data fft matrix of signals column by column
% Fe        : set of the effective harmonic numbers
% eff       : effective value of the signal X
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%

X = X(:);
X(1,:) = sqrt(2)*X(1,:);            % DC-component Correction
if nargin==2, X = X(Fe(:),:); end   % Input assignment
eff = sqrt(sum(abs(X).^2)/2);       % effective value

