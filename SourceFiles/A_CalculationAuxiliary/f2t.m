function x = f2t(X,N)
% F2T - Compute time domain signals from Fourier coefficients.
%
%    X    : Matrix containing the Fourier coefficients of
%           the different signals column by column
%    N    : Number of required time samples
%    x    : Matrix containing the time signals
% 
% See also TIME2FOUR.

if nargin==1
    [nRow,nCol] = size(X);
    Ndummy = nRow*2;
    N = 2;
    while Ndummy>N, N = N*2; end
end
x = N*real(ifft(X,N));

