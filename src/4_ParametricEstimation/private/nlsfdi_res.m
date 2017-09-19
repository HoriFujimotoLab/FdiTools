function cost = nlsfdi_res(Bn,An,freq,FRF,FRF_W,cORd,fs)
% NLSFDI_RES - Non Linear Least Squares FDI Residuals (MIMO).
%   cost = nlsfdi_res(Bn,An,freq,FRF,FRF_W,cORd,fs)
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
nroff = length(freq);           % number of frequency lines
nrofh = size(Bn,1);             % number of transfer functions
n=size(An,2)-1;                 % tranfer function order
nrofb = size(Bn,2);             % number of nominator parameters 

% calculation of frequency axis
j=sqrt(-1);
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   fprintf(' \n time domain undefined; it is set to continuous time \n')
   cORd = 'c';
   waxis = j*2*pi*freq;
end

% calculation of cost
EX = kron(ones(nroff,1),(n:-1:0));
W = kron(ones(1,n+1),waxis);
P = (W.^EX);

EX = kron(ones(nroff,1),(n:-1:-nrofb+n+1));
W = kron(ones(1,nrofb),waxis);
Q = (W.^EX);

E=[];
for h=1:nrofh
    E = [E; (FRF(:,h)-(Q*Bn(h,:)')./(P*An')).*FRF_W(:,h) ];
end
cost = (norm(E).^2)/2;

end