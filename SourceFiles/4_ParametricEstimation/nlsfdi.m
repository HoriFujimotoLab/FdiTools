function [Bn,An,Bwls,Awls]=nlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,iterno,relvar,GN,cORd,fs)
% NLLSFDI - Non-linear Least Squares FDI (MIMO).
%   [Bn,An,Bwls,Awls]=nlsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,iterno,relvar,GN,cORd,fs)
% FRF       : matrix of measured FRF values
% freq      : discete frequency vector
% FRF_W     : matrix of frequency weighting functions
% n         : order of the denominator polynomial
% M_mh,M_ml : vector of high and low order of the numerator polynomials
%             the i-th element of M_mh and M_ml corresponds to the i-th FRF
% iterno    : number of iterations
% relvar    : minimum relative deviation of the cost function
% Bn,An     : Non-linear least square solution
% Als,Bls   : Least square solution
% Awls,Bwls : Weighted least square solution
% GN        : GN=1 : Gauss Newton, GN=0: Levenberg-Marquardt
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model          
% fs        : sampling frequency (optional parameter)
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
% make vectors of M_mh and M_ml if they entered as matrices
M_mh=M_mh';
M_ml=M_ml';
M_mh = M_mh(:);
M_ml = M_ml(:);
nino=length(M_ml);
N = length(freq(:));                % number of FRF measurement points
nrofB=sum(M_mh-M_ml)+nino;          % total of numerator coefficients
nrofpar = nrofB+n;

% calculation of the frequency axis
j=sqrt(-1);
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   fprintf(' \n time domain is undefined; it is set to continuous time \n')
   cORd = 'c';
   waxis = j*2*pi*freq;
end;

if (max(M_mh) > n)
   fprintf(' \n The order of the numerator is larger than the order of the denominator \n');
   %return;
end;
if (min(M_mh-M_ml) < 0) 
   fprintf(' \n The elements of M_lh must be smaller that the corresponding elements of M_mh \n');
   return;
end;
   
% calculation of the frequency matrices used to avoid the use of polyval
ex=(n:-1:0)';
EX=kron(ones(N,1),ex');
W = kron(ones(1,n+1),waxis);
P_fr= (W.^EX);

Ntot = max(n+1,max(n-min(M_ml))+1);
ex=(n:-1:-Ntot+n+1)';
EX=kron(ones(N,1),ex');
W = kron(ones(1,Ntot),waxis);
Q_fr = (W.^EX);

%WLSFDI: initial estimate for iterative process
fprintf(' \n calculation of the LS solution \n')

FRF_WD = FRF.*FRF_W;
ex=(n:-1:0)';
EX=kron(ones(nino*N,1),ex');
W = kron(ones(nino,n+1),waxis);
P= (W.^EX).*kron(ones(1,n+1),FRF_WD(:));

Q = zeros(nino*N,nrofB);
index_count=1;
for (i=1:nino)
  U=kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_W(:,i)).*(kron(ones(1,M_mh(i)-M_ml(i)+1),waxis).^kron(ones(N,1),[M_mh(i):-1:M_ml(i)]));
  Q(N*(i-1)+1:N*i,index_count:index_count+M_mh(i)-M_ml(i))=U;
  index_count = index_count + M_mh(i)-M_ml(i)+1;
end; 

A = [real(P(:,2:n+1)) -real(Q);imag(P(:,2:n+1)) -imag(Q)];
b = -1*[real(P(:,1)) ; imag(P(:,1))];
y=pinv(A)*b;

% storing the solution in matrices An and Bn
% An is just a row vector
% Bn is a matrix, in which the polynomial coefficients are stored 
[Bn,An] = BA_construct(y,n,M_mh,M_ml);
Bls=Bn; Als=An;

% repeat the least squares estimate with better weighting functions
Ajw = P_fr*An';
FRF_WW = FRF_W./kron(ones(1,nino),abs(Ajw));
FRF_WD = FRF.*FRF_WW;

%LSFDI: initial estimate for iterative process
fprintf(' \n re-calculation of the LLS solution \n')

P= (W.^EX).*kron(ones(1,n+1),FRF_WD(:));
index_count=1;
for (i=1:nino)
  U=kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_WW(:,i)).*(kron(ones(1,M_mh(i)-M_ml(i)+1),waxis).^kron(ones(N,1),[M_mh(i):-1:M_ml(i)]));
  Q(N*(i-1)+1:N*i,index_count:index_count+M_mh(i)-M_ml(i))=U;
  index_count = index_count + M_mh(i)-M_ml(i)+1;
end
A = [real(P(:,2:n+1)) -real(Q);imag(P(:,2:n+1)) -imag(Q)];
b = -1*[real(P(:,1)) ; imag(P(:,1))];
y=pinv(A)*b;

% storing the solution in matrices An and Bn
% An is just a row vector
% Bn is a matrix, in which the polynomial coefficients are stored
[Bn,An] = BA_construct(y,n,M_mh,M_ml);
Awls = An; Bwls = Bn;

%calculation of cost for Bls, Als
cost0 = nlsfdi_res(Bwls,Awls,freq,FRF,FRF_W,cORd,fs);

% iterative estimation according to the Levenberg-Marquardt algorithm
fprintf('\n start of iteration \n');

if GN==1                    % initialization of parameters
  relax = 0;
else
 relax = 0.01;
end;

iter = 0;
iter0 = 0;
relerror0 = Inf;
relerror = Inf;
y0 = y;
dE_dB = zeros(nino*N,nrofB);

% start of iterative process
while (iter<iterno)&&(abs(relerror)>relvar)

  iter = iter + 1;
  Ajw = P_fr*An';               % Ajw = polyval(An,waxis);

  E = [];
  Bjw = [];
  for (i=1:nino)
    Bjw = [Bjw Q_fr*Bn(i,:)']; % Bjw = [Bjw polyval(Bn(i,:),waxis)];
    E = [E; (FRF(:,i)-Bjw(:,i)./Ajw).*FRF_W(:,i) ];
  end;

% scaling of the measured FRF with the weighting function
  FRF_WD = FRF_W./kron(ones(1,nino),Ajw.^2);

  dE_dA=[];
  if (n~=0)
    ex=(n-1:-1:0)';
    EX=kron(ones(nino*N,1),ex');
    W = kron(ones(nino,n),waxis);
    dE_dA= (W.^EX).*kron(ones(1,n),FRF_WD(:).*Bjw(:));
  end
  FRF_WD = FRF_W./kron(ones(1,nino),Ajw);

  index_count=1;
  for (i=1:nino)
    U=-kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_WD(:,i)).*(kron(ones(1,M_mh(i)-M_ml(i)+1),waxis).^kron(ones(N,1),[M_mh(i):-1:M_ml(i)]));
    dE_dB(N*(i-1)+1:N*i,index_count:index_count+M_mh(i)-M_ml(i))=U;
    index_count = index_count + M_mh(i)-M_ml(i)+1;
  end; 

  A = [real(dE_dA) real(dE_dB);imag(dE_dA) imag(dE_dB)];
  b = [real(E) ; imag(E)];
  JtJ = A'*A;
  Jte = A'*b;
  diagJtJ = diag(JtJ);
   
  A1 = A;
  A2 =  sqrt(relax*diag(diagJtJ+max(diagJtJ)*eps));
  A = [A1 ; A2];
  b = [b; zeros(nrofpar,1)];
  dy = - pinv(A)*b;  
  y0 = y;
  y = y + dy;

% storing the solution in matrices An and Bn
  [Bn,An]=BA_construct(y,n,M_mh,M_ml);
    
% calculation of the cost   
  cost = nlsfdi_res(Bn,An,freq,FRF,FRF_W,cORd,fs);
 
  relerror = (cost-cost0)/cost0; %relative deviation of the cost
  
  if ((cost < cost0)||(GN==1))
     y0 = y; cost0=cost; %updating solution 
     iter0 = iter;relerror0 = relerror;
     relax = relax/2; %lowering the Levenberg-Marquardt factor
  else
     y = y0;cost=cost0; %restoring the best result 
     [Bn,An] = BA_construct(y,n,M_mh,M_ml);
     relax = relax*10;
  end;    

  fprintf('teller = %g  ITER = %g  cost = %g  rel. change of cost = %g \n',iter,iter0,cost0,relerror0)
 
end;