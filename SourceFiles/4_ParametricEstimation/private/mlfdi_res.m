function cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,cORd,fs)
% MLFDI_RES - Maximum Likelihood Estimation Residuals.
%
% Bn,An     : solution after "iter" iterations
% freq      : discrete frequency vector
% X,Y       : input,output values of the FRF
% sX2       : variance of the input frequency domain noise 
% sY2       : variance of the output frequency domain noise 
% cXY       : covariance of input/output frequency domain noise
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model          
% fs        : sampling frequency (optional parameter)
% Author    : Thomas Beauduin, KULeuven, 2014
%%%%%
j=sqrt(-1);
freq=freq(:);
N=length(freq);
n=size(An,2)-1;

% calculation of the frequency axis
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   disp('time domain is undefined; it is set to continuous time');
   cORd = 'c';
   waxis = j*2*pi*freq;
end;

% matrices P and Q are used to form the complex set of equations
P = kron(ones(1,n+1),waxis).^kron(ones(N,1),(n:-1:0));
Num = P*Bn'; Den = P*An';
SE=sqrt(sX2.*(abs(Num).^2)+sY2.*(abs(Den).^2)-2*real(cXY.*Den.*conj(Num)));
E = (Num.*X - Den.*Y)./SE;
cost = (norm(E)^2)/2;

end
