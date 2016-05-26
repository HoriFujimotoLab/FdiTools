function X=schroed(X)
%SCHROED - Schroeder multisine phase design.
%
% X (in)    : abs(X) ampl. spectrum of the signals column by column.
% X (out)   : Matrix containing the spectrum with shroeder phases
% Algorithm : phase.n = phase.1 - 2p.[S{k=1,n-1} (n-k)Ak^2].
% Author    : Thomas Beauduin, KULeuven, 2014
%%%%%

ampl = abs(X);
[freqno,signo] = size(ampl);
phase = zeros(size(ampl));
amplnorm = ampl./(ones(size(freqno,1))*sqrt(sum(ampl.^2)));
amplnorm = 2*pi*amplnorm.^2;
for i=3:freqno,
   phase(i,:)=phase(i-1,:)-sum(amplnorm(1:i-1,:));
end
X=ampl.*exp(1i*phase);

