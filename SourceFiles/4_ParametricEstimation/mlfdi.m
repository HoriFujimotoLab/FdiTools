function [Hml,Hls] = mlfdi(X,Y,freq,sX2,sY2,cXY,n,M_mh,M_ml,iterno,relvar,GN,cORd,fs)
% MLFDI - Maximum Likelihood Estimation (MIMO).
%   [Hml,Hls]=mlfdi(X,Y,freq,sX2,sY2,cXY,n,M_mh,M_ml,iterno,relvar,GN,cORd,fs)
% X,Y,freq  : Input & output frequency domain data
% sX2,sY2   : variance of X & Y frequency domain data
% cXY       : Covariance between X & Y frequency domain data
% n,mh,ml   : Order of the denominator/nominator polynomials
% iterno    : Maximum number of iterations (stop criterion)
% relvar    : Maximum relative deviation of the cost function
% cORd, fs  : Continuous 'c' or discrete 'd' model identification
% Bm/l,Am/l : ML/LS iterative & initial estimation solution
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%
% See also WLSFDI, NLSFDI, GTLSFDI, BTLSFDI
%%%%%
nrofi = size(X,2);                      % number of inputs
nrofo = size(Y,2);                      % number of outputs
nrofh = nrofi*nrofo;                    % number of transfer functions
nroff = length(freq(:));                % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;           % number of numerator coefficients
nrofp = nrofb+n;                        % number of estimated parameters
M_mh=M_mh'; M_ml=M_ml';                 % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

% Calculation of initial values for iterative process
fprintf(' \n Initial calculation: LS solution \n')
[Hls,waxis] = lsfdi(X,Y,freq,n,M_mh,M_ml,cORd,fs);

% Calculation of iterative parameter estimation
fprintf('\n Iterative calculation: ML solution \n');
if GN==1,   relax = 0;                  % gradient relaxation
else        relax = 1;
end
[Bml,Aml] = hm2ba(Hls);                 % starting values choice
iter0 = 0; iter = 0;                    % interation number
relerror0 = Inf; relerror = Inf;        % relative error
y = ba2theta(Bml,Aml,n,M_mh,M_ml);      % initial parameter vector
dA = zeros(nrofh*nroff,n);              % denominator change
dB = zeros(nrofh*nroff,nrofb);          % numerator change
cost0 = mlfdi_res(Bml,Aml,freq,X,Y,sX2,sY2,cXY,waxis);

EX = kron(ones(nroff,1),(n:-1:0));
W = kron(ones(1,n+1),waxis);
P = (W.^EX);

EX = kron(ones(nroff,1),(max(M_mh):-1:min(M_ml)));
W = kron(ones(1,max(M_mh)-min(M_ml)+1),waxis);
Q = (W.^EX);

while (iter<iterno)&&(relerror>relvar)
  iter = iter + 1;
  Num = P*Bml'; Den = P*Aml';
  
  % Cost function minimization p365
  index = 0;
  E = zeros(nrofh*nroff,1);
  SE = zeros(nrofh*nroff,1);
  for h=1:nrofh
      i = ceil(h/nrofo); o = h-(i-1)*nrofo;
      SE((h-1)*nroff+1:h*nroff) = ...
            sqrt( sX2(:,i).*(abs(Num(:,h)).^2) ...
                + sY2(:,o).*(abs(Den).^2) ...
                - 2*real(cXY(:,h).*Den.*conj(Num(:,h))) );
      E((h-1)*nroff+1:h*nroff) = Num(:,h).*X(:,i) - Den.*Y(:,h);
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

  A = [J ; sqrt(relax*diag(diag(JtJ) + max(diag(JtJ))*eps))];
  b = [e; zeros(nrofp,1)];
  dy = -pinv(A)*b;   
  y0 = y;
  y = y + dy;

  [Bml,Aml] = theta2ba(y,n,M_mh,M_ml); 
  cost = mlfdi_res(Bml,Aml,freq,X,Y,sX2,sY2,cXY,waxis);
  relerror = abs(cost-cost0)/cost0;

  if ((cost < cost0)||(GN==1))
      y0 = y; cost0 = cost;             % updating solution 
      iter0 = iter;
      relerror0 = relerror;
      relax = relax/2;                  % lowering Levenberg factor
  else
      y = y0; cost = cost0;             % restoring best result 
      [Bml,Aml] = theta2ba(y,n,M_mh,M_ml);
      relax = relax*10;                 % increasing Levenberg factor
  end
  fprintf('Iter %g: index = %g, cost = %g, rel.err = %g\n',...
           iter,iter0,cost0,relerror0)   
end
Hml = ba2hm(Bml,Aml,nrofi,nrofo);

end