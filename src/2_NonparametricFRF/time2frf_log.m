function [Xi,Yi,FRF,freq,coh] = time2frf_log(x,y,fs,fl,fh,df,window,nrofl)
%TIME2FRF_LOG - logarithmic estimation of frf.
%
% x, y     : input and output measured data vector
% fs       : measurement sampling frequency
% fl, fh   : lowest and highest frequencies of excitation
% df       : spectral line differences
% window   : window choice: 0=rect,1=hann,2=cheb,3=tria
% nrofl    : number of overlap samples (welch windowing)
% FRF      : FRF-matrix [H11 H12 H13 .. H1m H21 ... Hlm]
% freq     : measured frequency lines
% coh      : multiple coherence
% Author   : Thomas Beauduin, KULeuven, 2014
%%%%%
[~,m] = size(x); 
[nroft,l] = size(y);
nrofs = fs/df;
freq = (fl:fh)'*df;

% Window Calculation
switch window
    case 0, win = ones(nrofs,1);    % rectangular window
    case 1, win = hanning(nrofs);   % hanning window
    case 2, win = chebwin(nrofs);   % Chebyshev window
    case 3, win = bartlett(nrofs);  % triangular window
end
winX = kron(win(:),ones(1,m));
winY = kron(win(:),ones(1,l));

% Windowed FFT
index = (1:nrofs);
k = fix((nroft-nrofl)/(nrofs-nrofl));
X = zeros(fh-fl+1,m*k); Y = zeros(fh-fl+1,l*k);
for i=1:k
   Xw = fft(x(index,:).*winX);
   Yw = fft(y(index,:).*winY);
   X(:,1+(i-1)*m:i*m) = Xw(fl+1:fh+1,:);
   Y(:,1+(i-1)*l:i*l) = Yw(fl+1:fh+1,:);
   index = index + (nrofs - nrofl);
end

% Calculate FRF & COH
FRF = zeros((fh-fl+1),m*l);
coh = zeros(fh-fl+1,l);
Xi = zeros(m,m); Yi=zeros(l,m);
for fk = 1:(fh-fl+1)
  Hi = zeros(l,m);
  for i=1:k-m+1
     for ii=1:m
       Xi(:,ii) = X(fk,1+(ii+i-2)*m:(ii+i-1)*m).'; 
       Yi(:,ii) = Y(fk,1+(ii+i-2)*l:(ii+i-1)*l).'; 
     end
     Hi =  Hi + Yi*inv(Xi);
  end
  Hi = unwrap(angle(Hi/(k-m+1)));
  Hav_fk = zeros(l,m);
  for i=1:k-m+1
    for ii=1:m
       Xi(:,ii) = X(fk,1+(ii+i-2)*m:(ii+i-1)*m).'; 
       Yi(:,ii) = Y(fk,1+(ii+i-2)*l:(ii+i-1)*l).'; 
    end
    Hav_fk = Hav_fk + log(Yi*inv(Xi).*exp(-1i*Hi));
  end
  Hav_fk = exp(Hav_fk/(k-m+1)).*exp(1i*Hi); 
  Hav_fk = Hav_fk.'; Hav_fk = Hav_fk(:);
  FRF(fk,:) = Hav_fk.';
  XiXi = Xi*Xi'; YiYi = Yi*Yi';
  for i=1:l
    coh(fk,i) = real((Hav_fk(i,:)*XiXi*Hav_fk(i,:)')/YiYi(i,i));
  end
end
end