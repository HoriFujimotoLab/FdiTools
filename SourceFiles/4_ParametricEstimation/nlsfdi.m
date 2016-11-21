function [Hnls,Hwls] = nlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,max_iter,max_err,GN,cORd,fs)
% NLLSFDI - Non-linear Least Squares FDI (MIMO).
%   [Hnls,Hwls] = nlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,max_iter,max_err,GN,cORd,fs)
% FRF,freq  : Transfer function frequency domain data
% FRF_W     : matrix of frequency weighting function
% n,mh,ml   : Order of the denominator/nominator polynomials
% max_iter  : Maximum number of iterations (stop criterion)
% max_err   : Maximum model relative error (stop criterion)
% GN        : GN=1 : Gauss Newton, GN=0: Levenberg-Marquardt
% cORd, fs  : Continuous 'c' or discrete 'd' model identification
% Bn/w,An/w : NLS/WLS iterative & initial estimation solution
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
nrofi = size(M_mh,2);           % number of inputs
nrofo = size(M_mh,1);           % number of outputs
nrofh = nrofi*nrofo;            % number of transfer functions
nroff = length(freq(:));        % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;   % number of numerator coefficients
nrofp = nrofb+n;                % number of estimated parameters

% Calculation of initial values for iterative process
fprintf(' \n Initial calculation: LS solution \n')
[Hwls,waxis] = wlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,cORd,fs);

M_mh = M_mh'; M_ml = M_ml';     % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

% Calculation of iterative parameter estimation
fprintf('\n Iterative calculation: NLS solution \n');
if GN==1,   relax = 0;                     % gradient relaxation
else        relax = 0.01;
end
dir = 1;
[Bn,An] = hm2ba(Hwls);                      % starting values choices
iter0 = 0; iter = 0;                        % interation number
relerror0 = Inf; relerror = Inf;            % relative error
dB = zeros(nrofh*nroff,nrofb);              % parameter error change
y = ba2theta(Bn,An,n,M_mh,M_ml);            % initial parameters
EX = kron(ones(nroff,1),(n:-1:0));          % External matrix
W = kron(ones(1,n+1),waxis);                % Weighting matrix
P_fr = (W.^EX); Q_fr = (W.^EX);             % First P/Q-matrix
cost0 = nlsfdi_res(Bn,An,freq,FRF,FRF_W,cORd,fs);

while (iter<max_iter)&&(abs(relerror)>max_err)
  iter = iter + 1;
  Ajw = P_fr*An';

  E = []; Bjw = [];
  for i=1:nrofh
      Bjw = [Bjw Q_fr*Bn(i,:)'];
      E = [E; (FRF(:,i)-Bjw(:,i)./Ajw).*FRF_W(:,i) ];
  end
  FRF_WD = FRF_W./kron(ones(1,nrofh),Ajw.^2);
  EX = kron(ones(nrofh*nroff,1),(n-1:-1:0));
  W = kron(ones(nrofh,n),waxis);
  dA = (W.^EX).*kron(ones(1,n),FRF_WD(:).*Bjw(:));
  
  index = 1;
  FRF_WD = FRF_W./kron(ones(1,nrofh),Ajw);
  for i=1:nrofh
      W = kron(ones(1,M_mh(i)-M_ml(i)+1),waxis);
      EX = kron(ones(nroff,1),[M_mh(i):-1:M_ml(i)]);
      U = (W.^EX).*-kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_WD(:,i));
      dB(nroff*(i-1)+1:nroff*i,index:index+M_mh(i)-M_ml(i)) = U;
      index = index + M_mh(i)-M_ml(i)+1;
  end
  J = [real(dA) real(dB) ; imag(dA) imag(dB)];
  e = [real(E) ; imag(E)];
  JtJ = J'*J;
  
  A = [J ; sqrt(relax*diag(diag(JtJ)+max(diag(JtJ))*eps))];
  b = [e; zeros(nrofp,1)];
  dy = -pinv(A)*b;  
  y0 = y;
  y = y + dy*dir;

  [Bn,An] = theta2ba(y,n,M_mh,M_ml); 
  cost = nlsfdi_res(Bn,An,freq,FRF,FRF_W,cORd,fs);
  relerror = (cost-cost0)/cost0;
  
  if ((cost < cost0)||(GN==1))
     y0 = y; cost0 = cost;          % updating solution 
     iter0 = iter; 
     relerror0 = relerror;
     relax = relax/2;               % lowering Levenberg factor
  else
     y = y0; cost = cost0;          % restoring best result 
     [Bn,An] = theta2ba(y,n,M_mh,M_ml);
     relax = relax*10;              % increasing Levenberg factor
  end
  if relax > 1e5
      dir = -1*dir;
      relax = 0.1;
  end
  
  fprintf('Iter %g: cost = %g, err = %g, r = %g \n',...
           iter, cost0, relerror0, relax)
end
Hnls = ba2hm(Bn,An,nrofi,nrofo);

end