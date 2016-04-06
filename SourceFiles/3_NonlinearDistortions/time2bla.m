function [X,Y,FRF,freq,Gbla,sX2,sY2,cXY] = time2bla(x,y,fs,fl,fh,df)
%TIME2BLA - Best Linear Approximation of FRF.
%   [X,Y,FRF,freq,Gbla,sX2,sY2,cXY] = time2bla(x,y,fs,fl,fh,df)
% x, y      : data vectors of periodic broad-band measurement
% fs,df     : sampling frequency and frequency resolution
% fl,fh     : lowest & highest frequency of excitated band
% FRFm,FRFn : measurement & noise frequency response functions
% sX2,sY2   : variance of real & imaginary parts of X,Y noise
% cXY       : covariance between real & imaginary parts of X,Y noise 
% sCR       : cramer-rao variance on measurement FRF
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
nrofs = fs/df;                                      % number of samples
nl = ceil(fl/df); nh = floor(fh/df);                % frequency range
freq = double((nl:1:nh)'/(nrofs/fs));               % frequency lines
nroff=length(freq);                                 % number of lines
nrofp = double(floor(length(x)/nrofs));             % number of periods
nrofm = min(size(x));                               % number of measures
if size(x,1) < size(x,2), x=x'; y=y'; end

% FFT
Xs=zeros(nrofm,nrofp,nrofs); %(M,P,S)
Ys=zeros(nrofm,nrofp,nrofs);
Fs=zeros(nrofm,nrofp,nrofs);
for m=1:nrofm
    for p=1:nrofp
        Xs(m,p,:)=fft(x(1+(p-1)*nrofs:p*nrofs,m));  % fft of 1x period
        Ys(m,p,:)=fft(y(1+(p-1)*nrofs:p*nrofs,m));
        Fs(m,p,:)=Ys(m,p,:)./Xs(m,p,:);
    end
end
Xf=Xs(:,:,nl+1:nh+1); %(M,P,F)                      % fft dc-term removal
Yf=Ys(:,:,nl+1:nh+1); 
Ff=Fs(:,:,nl+1:nh+1);

% BLA
Gbla(1,:) = squeeze(mean(mean(Ff,2),1));                             % mean
Gbla(2,:) = squeeze(std(mean(Ff, 2), 0, 1))/sqrt(nrofm);             % stdt
Gbla(3,:) = (squeeze(mean(std(Ff, 0, 2).^2, 1))/(nrofm*nrofp)).^0.5; % stdn
Gbla(4,:) = (nrofm*abs(Gbla(2,:).^2 - Gbla(3,:).^2)).^0.5;           % stds

% VAR
for m=1:nrofm
    sXp(m,:)=(squeeze(std(Xf(m,:,:)))'.^2)/2/nrofp; % measured variances
    sYp(m,:)=(squeeze(std(Yf(m,:,:)))'.^2)/2/nrofp;
    for i=1:nroff
        Q=cov(Xf(m,:,i),Yf(m,:,i));                 % measured covariance
        cXYp(m,i) = Q(1,2)/2/nrofp;
    end
end
sX2 = mean(sXp,1)'; %(M,F)->(F)
sY2 = mean(sYp,1)'; 
cXY = mean(cXYp,1)';

% FRF
X=squeeze(mean(mean(Xf,2),1)); 
Y=squeeze(mean(mean(Yf,2),1));                      % frequency averaging
FRF=Y./X;                                           % signal frf calc
sCR=2*abs(FRF).*(sX2./(abs(X)).^2 ...               % cramer-rao lowerbound
    +sY2./(abs(Y)).^2 ...
    -2*real(cXY./(conj(X).*Y)));

end
