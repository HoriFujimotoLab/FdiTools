function Fipp = cr_rao(X,Y,freq,B,A,sX2,sY2,cXY,n,mh,ml,cORd,fs)
%CR_RAO - Cramer-Rao Lower Bound of parameter covariance matrix.
%
% X         : input values of the FRF
% Y         : output values of the FRF
% freq      : frequency axis in Hz
% sX2       : variance of the input frequency domain noise 
% sY2       : variance of the output frequency domain noise 
% cXY       : covariance of input and output frequency domain noise 
% n         : order of the denominator polynomial
% mh,ml     : high and low order of the numerator polynomial
% B,A       : model in polynomial form 
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model
% fs        : sampling frequency (optional parameter)
% Fipp      : Fisher Information Matrix
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
B(1:n-mh) = [];
y = [A(2:n+1)';B(1:mh-ml+1)'];

j=sqrt(-1);
freq=freq(:);
N=length(freq);

% calculation of the frequency axis
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   fprintf('\n time domain is undefined; it is set to continuous time \n')
   cORd = 'c';
   waxis = j*2*pi*freq;
end;

P = kron(ones(1,n+1),waxis).^kron(ones(N,1),(n:-1:0));
Q = kron(ones(1,mh-ml+1),waxis).^kron(ones(N,1),(mh:-1:ml));

Num = Q*y(n+1:n+1+mh-ml);
Den = P*[1;y(1:n)];
SE = sqrt( sX2.*(abs(Num).^2) + sY2.*(abs(Den).^2) - 2*real(cXY.*Den.*conj(Num)) );
E = Num.*X - Den.*Y;

A = [];

for (i=2:n+1)
   WW = -Y.*P(:,i)./SE - E./(SE.^3).*( sY2.*real(Den.*conj(P(:,i))) ...
                                        - real(cXY.*P(:,i).*conj(Num)) );
   A = [A WW];
end;   


for (i=1:mh-ml+1)
   WW = X.*Q(:,i)./SE - E./(SE.^3).*( sX2.*real(Num.*conj(Q(:,i))) ...
                                       - real(cXY.*Den.*conj(Q(:,i))) );
   A = [A WW];
end;   

J=[real(A);imag(A)];

Fipp = J'*J;

end