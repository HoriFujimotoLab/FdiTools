function cost = btlsfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,relax,waxis)
% MLFDI_RES - Maximum Likelihood Estimation Residuals (MIMO).
%   cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,waxis)
% Bml,Aml   : ML iterative calculated estimation solution
% X,Y,freq  : Input & output frequency domain data
% sX2,sY2   : variance of X & Y frequency domain data
% cXY       : Covariance between X & Y frequency domain data
% relax     : weigthing relax-factor for 
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

E = zeros(nrofh*nroff,1);
SE2 = zeros(nrofh*nroff,1);
SEr2 = zeros(nrofh*nroff,1);
for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    SEr2((h-1)*nroff+1:h*nroff) = ...
            ( sX2(:,i).*(abs(Num(:,h)).^2)...
            + sY2(:,o).*(abs(Den).^2)...
            - 2*real(cXY(:,h).*Den.*conj(Num(:,h)))).^relax;
    SE2((h-1)*nroff+1:h*nroff) = ...
              sX2(:,i).*(abs(Num(:,h)).^2) ...
            + sY2(:,o).*(abs(Den).^2) ...
            - 2*real(cXY(:,h).*Den.*conj(Num(:,h)));
    E((h-1)*nroff+1:h*nroff) = abs(Num(:,h).*X(:,i) - Den.*Y(:,o)).^2;
end
cost = sum(E./SEr2)/sum(SE2./SEr2);

end