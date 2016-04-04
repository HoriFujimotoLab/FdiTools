function [Bg,Ag,X] = gtlsfdi(X,Y,freq,n,M_mh,M_ml,sX2,sY2,cXY,cORd,fs)
% GTLS - Generalized Total Least Squares Estimation (SISO).
%   [Bg,Ag,xqsvd] = gtlsfdi(Y,X,freq,n,mh,ml,sY2,sX2,cXY)
% X         : input values of the FRF
% Y         : output values of the FRF
% freq      : frequency vector
% sX2,sY2   : variance of input/output data
% cXY       : covariance of input/output data
% n         : order of the denominator polynomial
% mh,ml     : high and low order of the numerator polynomial
% Bg,Ag     : GTLS solution in polynomial form
% xqsvd     : QSVD solution (i.e. GTLS solution in vector form)
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
M_mh=M_mh'; M_ml=M_ml';             % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

nrofi = size(X,2);                  % number of inputs
nrofo = size(Y,2);                  % number of outputs
nrofh = nrofi*nrofo;                % number of transfer functions
nroff = length(freq(:));            % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;       % number of numerator coefficients
nrofp = nrofb+n;                    % number of estimated parameters

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

% Calculation of data matrix A
EX = kron(ones(nrofh*nroff,1),(n:-1:0));
W = kron(ones(nrofh,n+1),waxis);
for i=1:nrofi
    Yh(:,(i-1)*nrofo+1:i*nrofo) = Y;
end
P = (W.^EX).*kron(ones(1,n+1),Yh(:));

Q = zeros(nrofh*nroff,nrofb);
index = 1;
for h=1:nrofh
    EX = kron(ones(nroff,1),(M_mh(h):-1:M_ml(h)));
    W = kron(ones(1,M_mh(h)-M_ml(h)+1),waxis);
    U = (W.^EX).*kron(ones(1,M_mh(h)-M_ml(h)+1),X(:,ceil(h/nrofo)));
    Q(nroff*(h-1)+1:nroff*h,index:index+M_mh(h)-M_ml(h)) = U;
    index = index + M_mh(h)-M_ml(h)+1;
end
A = [real(P) -real(Q); imag(P) -imag(Q)];

% Calculation of right weighting matrix B
B = zeros(nrofh*(nrofp+1),nrofp+1);
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
    B(rows,cols) = chol(MtM);
end

% Generalized eigenvalue problem
X = qsvd(A,B);
X = inv(X');
X = X(:,nrofp+1)/X(1,nrofp+1);
[Bg,Ag] = BA_construct(X(2:end),n,M_mh,M_ml);

end