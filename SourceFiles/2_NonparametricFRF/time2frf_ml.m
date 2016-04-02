function [Xs,Ys,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,y,fs,fl,fh,df)
%TIME2FRF_ML - maximum-likelihood estimation of FRF (MIMO).
%   [Xs,Ys,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,y,fs,fl,fh,nrofs)
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
[~,nrofi] = size(x);                                % number of inputs
[~,nrofo] = size(y);                                % number of outputs 
nrofs = fs/df;                                      % samples per period
nl = ceil(fl/df); nh = floor(fh/df);                % low & high freq
freq = double((nl:1:nh)'/(nrofs/fs));               % full freq lines
nroff = length(freq);                               % number of freq lines
nrofp = double(floor(length(x)/nrofs));             % number of periods

% SIGNAL
INP = zeros(nroff,nrofp,nrofi); 
OUT = zeros(nroff,nrofp,nrofo);
Xs = zeros(nroff,nrofi); sX2 = zeros(nroff,nrofi);
Ys = zeros(nroff,nrofo); sY2 = zeros(nroff,nrofo);
cXY = zeros(nroff,nrofi,nrofo); 
FRFs = zeros(nroff,nrofi,nrofo); sCR = zeros(nroff,nrofi,nrofo);
for i=1:nrofi
    for p=1:nrofp
        Ip = fft(x(1+(p-1)*nrofs:p*nrofs,i));       % fft of 1x period
        INP(:,p,i) = Ip(nl+1:nh+1);                 % fft dc-term removal
    end
    Xs(:,i) = mean(INP(:,:,i),2);
    sX2(:,i)=((std(INP(:,:,i),0,2)).^2)/2/nrofp;    % measurement variances
end
for o=1:nrofo
    for p=1:nrofp
        Op = fft(y(1+(p-1)*nrofs:p*nrofs,o));
        OUT(:,p,o) = Op(nl+1:nh+1);
    end
    Ys(:,o) = mean(OUT(:,:,o),2);
    sY2(:,o)=((std(OUT(:,:,o),0,2)).^2)/2/nrofp;
end
for i=1:nrofi
    for o=1:nrofo
        for f=1:nroff
            Cf = cov(INP(f,:,i),OUT(f,:,o));        % measurement covariance
            cXY(f,i,o) = Cf(1,2)/2/nrofp;
        end
        FRFs(:,i,o) = Ys(:,o)./Xs(:,i);
        sCR(:,i,o) = 2*abs(FRFs(:,i,o)).*(sX2(:,i)./(abs(Xs(:,i))).^2 ...
            + sY2(:,o)./(abs(Ys(:,o))).^2 ...
            - 2*real(cXY(:,i,o)./(conj(Xs(:,i)).*Ys(:,o))));
    end
end

% NOISE
OUT = zeros(nroff*2,floor(nrofp/2),nrofo);
NSE = zeros(nroff,floor(nrofp/2),nrofo);
Yn = zeros(nroff,nrofo);
FRFn = zeros(nroff,nrofi,nrofo);
for o=1:nrofo
    for p=1:floor(nrofp/2)
        Op = fft(y(1+(p-1)*nrofs*2:p*nrofs*2,o));   % fft of 2x period
        OUT(:,p,o) = Op(2*nl:2*nh+1);               % fft dc-term removal
    end
    for f=2:2:nroff*2                               % uneven freq lines
        NSE(f/2,:,o) = OUT(f,:,o);
    end
    Yn(:,o) = mean(NSE(:,:,o),2);
end
for i=1:nrofi
    for o=1:nrofo
        FRFn(:,i,o) = Yn(:,o)./Xs(:,i);             % noise frf calc
    end
end
                                  
end