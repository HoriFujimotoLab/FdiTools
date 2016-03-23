function [Bg,Ag,xqsvd] = gtlsfdi(Y,X,freq,n,mh,ml,sY2,sX2,cXY)
% GTLS - Generalized Total Least Squares Estimation
%
% X         : input values of the FRF
% Y         : output values of the FRF
% freq      : frequency vector
% sX2       : variance of the input frequency domain noise 
% sY2       : variance of the output frequency domain noise 
% cXY       : covariance of between noise on X and Y
% n         : order of the denominator polynomial
% mh,ml     : high and low order of the numerator polynomial
% Bg,Ag     : GTLS solution in polynomial form
% xqsvd     : QSVD solution (i.e. GTLS solution in vector form)
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
freq=freq(:);
w = 1i*freq*2*pi;
nA = (n:-1:0);
nB = (mh:-1:ml);

A = (kron(w,ones(size(nA))).^kron(ones(size(w)),nA)).*kron(Y,ones(size(nA)));
A = [A -(kron(w,ones(size(nB))).^kron(ones(size(w)),nB)).*kron(X,ones(size(nB)))];
A = [real(A);imag(A)];

MytMy = zeros(n+1,n+1);
MxtMx = zeros(mh-ml+1,mh-ml+1);
MytMx = zeros(n+1,mh-ml+1);

for p=n:-1:0
  for q=n:-1:0
      MytMy(n-p+1,n-q+1) = real(((-1)^q)*2*sum(w.^(p+q).*sY2));
  end;
end;

for p=mh:-1:ml
  for q=mh:-1:ml
      MxtMx(mh-p+1,mh-q+1) = real(((-1)^q)*2*sum(w.^(p+q).*sX2));
  end;
end;

for p=n:-1:0
  for q=mh:-1:ml
      MytMx(n-p+1,mh-q+1) =  -real(((-1)^q)*2*sum(w.^(p+q).*cXY));
  end;
end;

MtM=[MytMy MytMx ; MytMx' MxtMx];
M = chol(MtM);
  
%qsvd solution
[Ua,Un,Sa,Sn,xqsvd]=qsvd(A,M);
nt=n+1+mh-ml+1;
xqsvd=inv(xqsvd');
xqsvd=xqsvd(:,nt)/xqsvd(1,nt);

Ag = xqsvd(1:n+1)';
Bg = [zeros(1,n-mh) xqsvd(n+2:nt)' zeros(1,ml)];

end
