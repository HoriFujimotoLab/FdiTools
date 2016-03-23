function [Bls,Als]=lsfdi(FRF,freq,FRF_W,n,M_mh,M_ml,cORd,fs)
%LLSFDI - Linear Least Squares FDI (MIMO)
%
% FRF       : matrix of measured FRF values
% freq      : measured frequency lines vector
% FRF_W     : matrix of frequency weighting function
% n         : order of the denominator polynomial
% M_mh,M_ml : high and low order of the numerator polynomials
% Bls,Als   : LS-solution
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model         
% fs        : sampling frequency (optional parameter)
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
M_mh=M_mh'; M_ml=M_ml';
M_mh = M_mh(:); M_ml = M_ml(:);

nino=length(M_ml);
N = length(freq(:));
nrofB=sum(M_mh-M_ml)+nino;  % total of numerator coefficients
nrofpar = nrofB+n; 
j=sqrt(-1);

% calculation of the frequency axis
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   fprintf(' \n time domain is undefined; it is set to continuous time \n')
   cORd = 'c';
   waxis = j*2*pi*freq;
end

% scaling of the measured FRF with the weighting function
FRF_WD = FRF.*FRF_W;

%LLS-estimation
ex=(n:-1:0)';
EX=kron(ones(nino*N,1),ex')
W = kron(ones(nino,n+1),waxis);
P= (W.^EX).*kron(ones(1,n+1),FRF_WD(:));
Q = zeros(nino*N,nrofB);
index_count=1;
for (i=1:nino)
  U=kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_W(:,i)).*(kron(ones(1,M_mh(i)-M_ml(i)+1),waxis).^kron(ones(N,1),[M_mh(i):-1:M_ml(i)]));
  Q(N*(i-1)+1:N*i,index_count:index_count+M_mh(i)-M_ml(i))=U;
  index_count = index_count + M_mh(i)-M_ml(i)+1;
end 
A = [real(P(:,2:n+1)) -real(Q);imag(P(:,2:n+1)) -imag(Q)];
b = -1*[real(P(:,1)) ; imag(P(:,1))];
y=pinv(A)*b;

% storing the solution in matrices An and Bn
[Bls,Als] = BA_const(y,n,M_mh,M_ml);

end