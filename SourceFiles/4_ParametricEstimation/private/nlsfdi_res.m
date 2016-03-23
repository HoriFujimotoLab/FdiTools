function cost = nlsfdi_res(Bn,An,freq,FRF,FRF_W,cORd,fs)
% NLSFDI_RES - Non Linear Least Squares FDI Residuals
%
% Bn,An     : Non-linear least square solution
% freq      : discete frequency vector
% FRF       : matrix of measured FRF values
% FRF_W     : matrix of frequency weighting functions
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model
% fs        : sampling frequency (optional parameter)
% cost      : residual cost function value
% Author    : Thomas Beauduin, KULeuven, 2014
%%%%%
j=sqrt(-1);
N = length(freq);
[nino,x] = size(Bn);

% calculation of the frequency axis
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   disp('time domain is undefined; it is set to continuous time');
   cORd = 'c';
   waxis = j*2*pi*freq;
end

% cost calculation
[x,n]=size(An);
n=n-1;
ex=(n:-1:0)';
EX=kron(ones(N,1),ex');
W=kron(ones(1,n+1),waxis);
P=(W.^EX);
Ajw = P*An';

[x,Ntot]=size(Bn);
ex=(n:-1:-Ntot+n+1)';
EX=kron(ones(N,1),ex');
W = kron(ones(1,Ntot),waxis);
Q = (W.^EX);

E=[];
for (i=1:nino)
  E = [E; (FRF(:,i)-(Q*Bn(i,:)')./Ajw).*FRF_W(:,i) ];
end
cost = sum(abs(E).^2)/2;

end