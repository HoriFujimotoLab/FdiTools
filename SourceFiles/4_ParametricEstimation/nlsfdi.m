function [Bn,An,Bwls,Awls] = nlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,max_iter,max_err,GN,cORd,fs)
% NLLSFDI - Non-linear Least Squares FDI (MIMO).
%   [Bn,An,Bwls,Awls]=nlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,max_iter,max_err,GN,cORd,fs)
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
M_mh=M_mh'; M_ml=M_ml';         % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

nrofh = length(M_ml);           % number of transfer functions (2)
nroff = length(freq(:));        % number of frequency lines (240)
nrofb = sum(M_mh-M_ml)+nrofh;   % number of numerator coefficients (4)
nrofp = nrofb+n;                % number of estimated parameters (8)

% INITIAL WLSFDI 
% initial values estimate for iterative process
fprintf(' \n Initial calculation: LS solution \n')
[~,Als,waxis] = wlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,cORd,fs);

EX = kron(ones(nroff,1),(n:-1:0));
W = kron(ones(1,n+1),waxis);
P_fr = (W.^EX);                             % First P-matrix
Q_fr = (W.^EX);                             % First Q-matrix

Ajw = P_fr*Als';
FRF_WW = FRF_W./kron(ones(1,nrofh),abs(Ajw)); % improved weigthing
[Bwls,Awls,waxis] = wlsfdi(FRF,freq,FRF_WW,n,M_mh,M_ml,cORd,fs);

% ITERATIVE ESTIMATION
% Levenberg-Marquardt algorithm NLS estimation
fprintf('\n Iterative calculation: NLS solution \n');
if GN==1,   relax = 0;                      % gradient relaxation
else        relax = 0.01;
end
An = Awls; Bn = Bwls;                       % choice of starting values
iter0 = 0; iter = 0;                        % interation number
relerror0 = Inf; relerror = Inf;            % relative error
dE_dB = zeros(nrofh*nroff,nrofb);           % parameter error change
y = ba2theta(Bn,An,n,M_mh,M_ml);             % initial parameters
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

  dE_dA = [];
  EX = kron(ones(nrofh*nroff,1),(n-1:-1:0));
  W = kron(ones(nrofh,n),waxis);
  dE_dA = (W.^EX).*kron(ones(1,n),FRF_WD(:).*Bjw(:));
  FRF_WD = FRF_W./kron(ones(1,nrofh),Ajw);

  index = 1;
  for i=1:nrofh
      W = kron(ones(1,M_mh(i)-M_ml(i)+1),waxis);
      EX = kron(ones(nroff,1),[M_mh(i):-1:M_ml(i)]);
      U = (W.^EX).*-kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_WD(:,i));
      dE_dB(nroff*(i-1)+1:nroff*i,index:index+M_mh(i)-M_ml(i)) = U;
      index = index + M_mh(i)-M_ml(i)+1;
  end
  J = [real(dE_dA) real(dE_dB) ; imag(dE_dA) imag(dE_dB)];
  e = [real(E) ; imag(E)];
  JtJ = J'*J;
  Jte = J'*e;
  diagJtJ = diag(JtJ);
  
  A1 = J;
  A2 = sqrt(relax*diag(diagJtJ+max(diagJtJ)*eps));
  A = [A1 ; A2];
  b = [e; zeros(nrofp,1)];
  dy = - pinv(A)*b;  
  y0 = y;
  y = y + dy;

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
  fprintf('Iter %g: index = %g, cost = %g, rel.err = %g\n',...
           iter, iter0, cost0, relerror0)
 
end