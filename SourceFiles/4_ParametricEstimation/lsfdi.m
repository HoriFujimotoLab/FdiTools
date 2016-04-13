function [Hls,waxis] = lsfdi(X,Y,freq,n,M_mh,M_ml,cORd,fs)
%LLSFDI - Linear Least Squares FDI (MIMO).
%   [Hls,waxis] = lsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,cORd,fs)
% X,Y,freq  : Input & output frequency domain data
% n,mh,ml   : Order of the denominator/nominator polynomials
% cORd, fs  : Continuous 'c' or discrete 'd' model identification
% Bm/l,Am/l : ML/LS iterative & initial estimation solution
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
nrofi = size(X,2);                  % number of inputs
nrofo = size(Y,2);                  % number of outputs
nrofh = nrofi*nrofo;                % number of transfer functions
nroff = length(freq(:));            % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;       % number of numerator coefficients
M_mh=M_mh'; M_ml=M_ml';             % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

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

% Calculation of jacobian matrix J
EX = kron(ones(nrofh*nroff,1),(n:-1:0));
W = kron(ones(nrofh,n+1),waxis);
for i=1:nrofi
    Yh(:,(i-1)*nrofo+1:i*nrofo) = Y;
end
P = (W.^EX).*kron(ones(1,n+1),Yh(:));

Q = zeros(nrofh*nroff,nrofb);
index = 1;
for h=1:nrofh
    EX = kron(ones(nroff,1),(M_mh(h):-1:M_ml(h)));
    W = kron(ones(1,M_mh(h)-M_ml(h)+1),waxis);
    U = (W.^EX).*kron(ones(1,M_mh(h)-M_ml(h)+1),X(:,ceil(h/nrofo)));
    Q(nroff*(h-1)+1:nroff*h,index:index+M_mh(h)-M_ml(h)) = U;
    index = index + M_mh(h)-M_ml(h)+1;
end
J = [real(P(:,2:n+1)) -real(Q);imag(P(:,2:n+1)) -imag(Q)];

% Calculation of pseudoinverse least squares solution
b = -1*[real(P(:,1)) ; imag(P(:,1))];
y = pinv(J)*b;
[Bls,Als] = theta2ba(y,n,M_mh,M_ml);
Hls = ba2hm(Bls,Als,nrofi,nrofo);

end