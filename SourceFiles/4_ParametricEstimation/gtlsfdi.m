function [Hgtls,waxis] = gtlsfdi(X,Y,freq,n,M_mh,M_ml,sX2,sY2,cXY,cORd,fs)
% GTLS - Generalized Total Least Squares Estimation (MIMO).
%   [Hgtls,waxis] = gtlsfdi(X,Y,freq,n,M_mh,M_ml,sX2,sY2,cXY,cORd,fs)
% X,Y,freq  : Input & output frequency domain data
% sX/Y2,cXY : Measurement frequency covariance matrix data
% n,mh,ml   : Order of the denominator/nominator polynomials
% cORd, fs  : Continuous 'c' or discrete 'd' model identification
% Bm/l,Am/l : ML/LS iterative & initial estimation solution
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
nrofi = size(X,2);                  % number of inputs
nrofo = size(Y,2);                  % number of outputs
nrofh = nrofi*nrofo;                % number of transfer functions
nroff = length(freq(:));            % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;       % number of numerator coefficients
nrofp = nrofb+n;                    % number of estimated parameters
M_mh=M_mh'; M_ml=M_ml';             % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

% Calculation of frequency axis
if     (cORd == 'c'),   waxis = 1i*2*pi*freq;
elseif (cORd == 'd'),   waxis = exp(1i*2*pi*freq/fs);
else                    waxis = 1i*2*pi*freq;
    fprintf('\n Time domain undefined; it is set to continuous time \n');                       
end
if (max(M_mh) > n)
   fprintf('\n Warning: numerator order is larger than denominator order\n');
end
if (min(M_mh-M_ml) < 0) 
   fprintf('\n Error: elements of M_ml must be smaller that those of M_mh\n');
   return;
end

% Calculation of weighted Jacobian (WreJre(Z))
Xh = zeros(nroff,nrofh);
Yh = zeros(nroff,nrofh);
for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    Xh(:,h) = X(:,i); Yh(:,h) = Y(:,o);
end
EX = kron(ones(nrofh*nroff,1),(n:-1:0));
W = kron(ones(nrofh,n+1),waxis);
P = (W.^EX).*kron(ones(1,n+1),Yh(:));

index = 1;
Q = zeros(nrofh*nroff,nrofb);
for h=1:nrofh
    EX = kron(ones(nroff,1),(M_mh(h):-1:M_ml(h)));
    W = kron(ones(1,M_mh(h)-M_ml(h)+1),waxis);
    U = (W.^EX).*kron(ones(1,M_mh(h)-M_ml(h)+1),Xh(:,h));
    Q(nroff*(h-1)+1:nroff*h,index:index+M_mh(h)-M_ml(h)) = U;
    index = index + M_mh(h)-M_ml(h)+1;
end
J = [real(P) -real(Q); imag(P) -imag(Q)];

% Calculation of column covariance matrix (Cwj^1/2)
C = zeros(nrofh*(nrofp+1),nrofp+1);
for h=1:nrofh
    i=ceil(h/nrofo); o=h-(i-1)*nrofo;
    MytMy = zeros(n+1,n+1);
    MxtMx = zeros(M_mh(h)-M_ml(h)+1,M_mh(h)-M_ml(h)+1);
    MytMx = zeros(n+1,M_mh(h)-M_ml(h)+1);
    for p=n:-1:0
        for q=n:-1:0
            MytMy(n-p+1,n-q+1) = ...
            2*real(((-1)^q)*sum(waxis.^(p+q).*sY2(:,o)));
        end
    end
    for p=M_mh(h):-1:M_ml(h)
        for q=M_mh(h):-1:M_ml(h)
            MxtMx(M_mh(h)-p+1,M_mh(h)-q+1) = ...
            2*real(((-1)^q)*sum(waxis.^(p+q).*sX2(:,i)));
        end
    end
    for p=n:-1:0
        for q=M_mh(h):-1:M_ml(h)
            MytMx(n-p+1,M_mh(h)-q+1) = ...
            -2*real(((-1)^q)*sum(waxis.^(p+q).*cXY(:,h)));
        end
    end
    MtM=[MytMy MytMx ; MytMx' MxtMx];
    cols = (1:n+1+M_mh(h)-M_ml(h)+1);
    rows = (h-1)*(nrofp+1)+1:(h-1)*(nrofp+1)+(n+1+M_mh(h)-M_ml(h)+1);
    C(rows,cols) = chol(MtM);
end

% Calculation of generalized right singular vector (Xg)
Xg = qsvd(J,C);
Xg = inv(Xg');
Xg = Xg(:,nrofp+1)/Xg(1,nrofp+1);
[Bgtls,Agtls] = theta2ba(Xg(2:end),n,M_mh,M_ml);
Hgtls = ba2hm(Bgtls,Agtls,nrofi,nrofo);

end