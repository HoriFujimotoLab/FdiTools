function [Bwls,Awls,waxis] = wlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,cORd,fs)
%LLSFDI - Linear Least Squares FDI (MIMO).
%   [Bls,Als]=lsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,cORd,fs)
% FRF       : matrix of measured FRF values
% freq      : measured frequency lines vector
% FRF_W     : matrix of frequency weighting function
% n         : order of the denominator polynomial
% M_mh,M_ml : high and low order of the numerator polynomials
% Bls,Als   : LS-solution (A - row vector, B - matrix)
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model         
% fs        : sampling frequency (optional parameter)
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
M_mh=M_mh'; M_ml=M_ml';         % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

nrofh = length(M_ml);           % number of transfer functions
nroff = length(freq(:));        % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;   % number of numerator coefficients

% Calculation of frequency axis
if     (cORd == 'c'),   waxis = 1i*2*pi*freq;
elseif (cORd == 'd'),   waxis = exp(1i*2*pi*freq/fs);
else                    waxis = 1i*2*pi*freq;
    fprintf('\n Time domain undefined; it is set to continuous time \n');                       
end
if (max(M_mh) > n)
   fprintf('\n Warning: numerator order is larger than denominator order\n');
end
if (min(M_mh-M_ml) < 0) 
   fprintf('\n Error: elements of M_ml must be smaller that those of M_mh\n');
   return;
end

% Calculation of data matrix A
FRF_WD = FRF.*FRF_W;
EX = kron(ones(nrofh*nroff,1),(n:-1:0));
W = kron(ones(nrofh,n+1),waxis);
P = (W.^EX).*kron(ones(1,n+1),FRF_WD(:));

Q = zeros(nrofh*nroff,nrofb);
index = 1;
for i=1:nrofh
  EX = kron(ones(nroff,1),(M_mh(i):-1:M_ml(i)));
  W = kron(ones(1,M_mh(i)-M_ml(i)+1),waxis);
  U = (W.^EX).*kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_W(:,i));
  Q(nroff*(i-1)+1:nroff*i,index:index+M_mh(i)-M_ml(i)) = U;
  index = index + M_mh(i)-M_ml(i)+1;
end 
A = [real(P(:,2:n+1)) -real(Q);imag(P(:,2:n+1)) -imag(Q)];

% Calculation of pseudoinverse least squares solution
b = -1*[real(P(:,1)) ; imag(P(:,1))];
y = pinv(A)*b;
[Bwls,Awls] = BA_construct(y,n,M_mh,M_ml);

end