function cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,waxis)
% MLFDI_RES - Maximum Likelihood Estimation Residuals (MIMO).
%   cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,waxis)
% Bml,Aml   : ML iterative calculated estimation solution
% X,Y,freq  : Input & output frequency domain data
% sX2,sY2   : variance of X & Y frequency domain data
% cXY       : Covariance between X & Y frequency domain data
% waxis     : Continuous or discrete time frequency axis 
% cost      : Maximum likelihood estimation residual cost
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
n = size(An,2)-1;                   % system order
nrofi = size(X,2);                  % number of inputs
nrofo = size(Y,2);                  % number of outputs
nroff = length(freq);               % number of frequency lines
nrofh = nrofi*nrofo;                % number of transfer functions

% calculation of residual estimation cost
EX = kron(ones(nroff,1),(n:-1:0));
W = kron(ones(1,n+1),waxis);
P = (W.^EX);
Num = P*Bn'; Den = P*An';

E = [];
for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    SE = sqrt(sX2(:,i).*(abs(Num(:,h)).^2) + ...
              sY2(:,o).*(abs(Den).^2) - ...
              2*real(cXY(:,h).*Den.*conj(Num(:,h))));
    E = [E; (Num(:,h).*X(:,i) - Den.*Y(:,o))./SE];
end
cost = (norm(E).^2)/2;

end