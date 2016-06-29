function X=t2f(x,N)
%T2F Fourier serie coefficients from time domain signals.
%
% x     : Matrix containing the time signals
% N     : Number of required time samples
% X     : Matrix of Fourier coefficients
% Author: Thomas Beauduin, KULeuven 2014

if nargin==1
   [rowno,colno]=size(x);
   Ndummy=max(rowno,colno);
   N=2;
   while Ndummy>N, N=N*2; end
end
X = fft(x,N);
X = [X(1,:);2*X(2:floor(N/2),:)]/N;

