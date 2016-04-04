function [Bn,An,Bls,Als,cost0]=mlfdi(X,Y,freq,sX2,sY2,cXY,n,M_mh,M_ml,iterno,relvar,GN,cORd,fs)
% MLFDI - Maximum Likelihood Estimation (MIMO).
%   [Bn,An,Bls,Als,cost0]=mlfdi(X,Y,freq,sX2,sY2,cXY,n,M_mh,M_ml,iterno,relvar,GN,cORd,fs)
% X,Y       : input,output values of the FRF
% freq      : discrete frequency vector
% sX2       : variance of the input frequency domain noise 
% sY2       : variance of the output frequency domain noise 
% cXY       : covariance of input and output frequency domain noise 
% n         : order of the denominator polynomial
% mh,ml     : high and low order of the numerator polynomial
% iterno    : number of iterations
% relvar    : minimum relative deviation of the cost function
% Bn,An     : solution after "iter" iterations
% Bls,Als   : LS-solution
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model          
% fs        : sampling frequency (optional parameter)
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
M_mh=M_mh'; M_ml=M_ml';                 % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

nrofi = size(X,2);                      % number of inputs
nrofo = size(Y,2);                      % number of outputs
nrofh = nrofi*nrofo;                    % number of transfer functions
nroff = length(freq(:));                % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;           % number of numerator coefficients
nrofp = nrofb+n;                        % number of estimated parameters

% Calculation of initial values for iterative process
fprintf(' \n Initial calculation: LS solution \n')
[Bls,Als,waxis] = lsfdi(X,Y,freq,n,M_mh,M_ml,cORd,fs);

% Calculation of iterative parameter estimation
fprintf('\n Iterative calculation: ML solution \n');
if GN==1,   relax = 0;                  % gradient relaxation
else        relax = 1;
end
An = Als; Bn = Bls;                     % starting values choice
iter0 = 0; iter = 0;                    % interation number
relerror0 = Inf; relerror = Inf;        % relative error
y = ba2yvec(Bn,An,n,M_mh,M_ml);         % initial parameter vector
dA = zeros(nrofh*nroff,n);              % denominator change
dB = zeros(nrofh*nroff,nrofb);          % numerator change
cost0 = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,cORd,fs);

EX = kron(ones(nroff,1),(n:-1:0));
W = kron(ones(1,n+1),waxis);
P = (W.^EX);

EX = kron(ones(nroff,1),(max(M_mh):-1:min(M_ml)));
W = kron(ones(1,max(M_mh)-min(M_ml)+1),waxis);
Q = (W.^EX);

while (iter<iterno)&&(relerror>relvar)
  iter = iter + 1;
  Num = P*Bn'; Den = P*An';
  
  E = []; SE = []; index = 0;
  for h=1:nrofh
      i = ceil(h/nrofo); o = h-(i-1)*nrofo;
      SE = [SE ; sqrt( sX2(:,i).*(abs(Num(:,h)).^2) ...
                 + sY2(:,o).*(abs(Den).^2) ...
                 - 2*real(cXY(:,h).*Den.*conj(Num(:,h))) )];
      E = [E ; Num(:,h).*X(:,i) - Den.*Y(:,h)];
      for j=2:n+1
          WW = - Y(:,o).*P(:,j)./SE((h-1)*nroff+1:h*nroff) ...
               - E((h-1)*nroff+1:h*nroff)./(SE((h-1)*nroff+1:h*nroff).^3)...
                 .*(sY2(:,o).*real(Den.*conj(P(:,j)))...
               - real(cXY(:,h).*P(:,j).*conj(Num(:,h))));
          dA(nroff*(h-1)+1:nroff*h,j-1) = WW;
      end
      for j=1:M_mh(h)-M_ml(h)+1
          WW =   X(:,i).*Q(:,j)./SE((h-1)*nroff+1:h*nroff) ...
               - E((h-1)*nroff+1:h*nroff)./(SE((h-1)*nroff+1:h*nroff).^3)...
                 .*(sX2(:,i).*real(Num(:,h).*conj(Q(:,j)))...
               - real(cXY(:,h).*Den.*conj(Q(:,j))) );
          dB(nroff*(h-1)+1:nroff*h,j+index) = WW;
      end
      index = index + M_mh(h)-M_ml(h)+1;
  end
  J = [real(dA) real(dB) ; imag(dA) imag(dB)];
  e = [real(E)./SE ; imag(E)./SE];
  JtJ = J'*J;
  Jte = J'*e;
  diagJtJ = diag(JtJ);

  A1 = J;
  A2 = sqrt(relax*diag(diagJtJ+max(diagJtJ)*eps));
  A = [A1 ; A2];
  b = [e; zeros(nrofp,1)];
  dy = -pinv(A)*b;   
  y0 = y;
  y = y + dy;

  [Bn,An] = BA_construct(y,n,M_mh,M_ml); 
  cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,cORd,fs);
  relerror = abs(cost-cost0)/cost0;

  if ((cost < cost0)||(GN==1))
      y0 = y; cost0 = cost;             % updating solution 
      iter0 = iter;
      relerror0 = relerror;
      relax = relax/2;                  % lowering Levenberg factor
  else
      y = y0; cost = cost0;             % restoring best result 
      [Bn,An] = BA_construct(y,n,M_mh,M_ml);
      relax = relax*10;                 % increasing Levenberg factor
  end
  fprintf('Index = %g iter = %g cost = %g rel.error = %g \n',...
          iter,iter0,cost0,relerror0)
      
end