function cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,cORd,fs)
% MLFDI_RES - Maximum Likelihood Estimation Residuals.
%
% Bn,An     : solution after "iter" iterations
% freq      : discrete frequency vector
% X,Y       : input,output values of the FRF
% sX2       : variance of the input frequency domain noise 
% sY2       : variance of the output frequency domain noise 
% cXY       : covariance of input/output frequency domain noise
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model          
% fs        : sampling frequency (optional parameter)
% Author    : Thomas Beauduin, KULeuven, 2014
%%%%%
nroff = length(freq);           % number of frequency lines
nrofh = size(Bn,1);             % number of transfer functions
n = size(An,2)-1;               % tranfer function order

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
Num = P*Bn'; Den = P*An';

E = [];
for i=1:nrofh
    SE = sqrt(sX2(:,i).*(abs(Num(:,i)).^2) + ...
              sY2(:,i).*(abs(Den).^2) - ...
              2*real(cXY(:,i).*Den.*conj(Num(:,i))));
    E = [E; (Num(:,i).*X(:,i) - Den.*Y(:,i))./SE];
end
cost = (norm(E).^2)/2;

end