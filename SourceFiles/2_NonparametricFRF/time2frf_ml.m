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
% Author    : Thomas Beauduin, KULeuven, PMA, Feb-2014
%%%%%
% future release will include full MIMO FRF processing
% future release will include structured input/output:
% EXC (data)
% <>.fs         :
% <>.fl         :
% <>.fh         :
% <>.df         :
% <>.ex         :
% TIME DATA
% <>.x          :
% <>.y          :
% <>.fl         :
% <>.fh         :
% FREQ DATA
% <>.freq       :
% <>.FRFs       :
% <>.X          :
% <>.Y          :
% <>.fs         : [Hz]
% COV
% <>.sX2        :
% <>.sY2        :
% <>.cXY        :
% <>.FRFn       :
% <>.sCR        :
% OPTIONS
% <>.frf_type   : ml/lpm/lrm/log/h2/h1
%%%%

[~,nrofi] = size(x);                                % number of inputs
[~,nrofo] = size(y);                                % number of outputs
nrofh = nrofi*nrofo;                                % number of tf's (Hxy)
nrofs = fs/df;                                      % samples per period
nl = ceil(fl/df); nh = floor(fh/df);                % low & high freq
freq = double((nl:1:nh)'/(nrofs/fs));               % full freq lines
nroff = length(freq);                               % number of freq lines
nrofp = double(floor(length(x)/nrofs));             % number of period

% Calculation of signal fft data
INP = zeros(nroff,nrofp,nrofi); 
OUT = zeros(nroff,nrofp,nrofo);
Xs = zeros(nroff,nrofi); sX2 = zeros(nroff,nrofi);
Ys = zeros(nroff,nrofo); sY2 = zeros(nroff,nrofo);
cXY = zeros(nroff,nrofh); 
FRFs = zeros(nroff,nrofh); 
sCR = zeros(nroff,nrofh);
for i=1:nrofi
    for p=1:nrofp
        Ip = fft(x(1+(p-1)*nrofs:p*nrofs,i));       % fft of 1x period
        INP(:,p,i) = Ip(nl+1:nh+1);                 % fft dc-term removal
 %       for f=1:nroff
 %           INP(f,p,i)=INP(f,p,i)*( 1i*2*pi*freq(f)/(1-exp(1i*2*pi*freq(f)/fs)) );
 %       end
    end
    Xs(:,i) = mean(INP(:,:,i),2);
    sX2(:,i)=((std(INP(:,:,i),0,2)).^2)/2/nrofp;    % measurement variances
end
for o=1:nrofo
    for p=1:nrofp
        Op = fft(y(1+(p-1)*nrofs:p*nrofs,o));
        OUT(:,p,o) = Op(nl+1:nh+1);
  %      for f=1:nroff
  %          INP(f,p,o)=INP(f,p,o)*( 1i*2*pi*freq(f)/(1-exp(1i*2*pi*freq(f)/fs)) );
  %      end
    end
    Ys(:,o) = mean(OUT(:,:,o),2);
    sY2(:,o)=((std(OUT(:,:,o),0,2)).^2)/2/nrofp;
end
for i=1:nrofi
    for o=1:nrofo
        for f=1:nroff
            Cf = cov(INP(f,:,i),OUT(f,:,o));        % measurement covariance
            cXY(f,(i-1)*nrofo+o) = Cf(1,2)/2/nrofp;
        end
        FRFs(:,(i-1)*nrofo+o) = Ys(:,o)./Xs(:,i);
        sCR(:,(i-1)*nrofo+o) = ...
            2*abs(FRFs(:,(i-1)*nrofo+o)).*(sX2(:,i)./(abs(Xs(:,i))).^2 ...
            + sY2(:,o)./(abs(Ys(:,o))).^2 ...
            - 2*real(cXY(:,(i-1)*nrofo+o)./(conj(Xs(:,i)).*Ys(:,o))));
    end
end

% Calculation of noise fft data
OUT = zeros(nroff*2,floor(nrofp/2),nrofo);
NSE = zeros(nroff,floor(nrofp/2),nrofo);
Yn = zeros(nroff,nrofo); FRFn = zeros(nroff,nrofh);
for o=1:nrofo
    for p=1:floor(nrofp/2)
        Op = fft(y(1+(p-1)*nrofs*2:p*nrofs*2,o));   % fft of 2x period
        OUT(:,p,o) = Op(2*nl:2*nh+1);               % fft dc-term removal
    end
    index = 1;
    for f=1:2:nroff*2                               % uneven freq lines
        NSE(index,:,o) = OUT(f,:,o);
        index = index + 1;
    end
    Yn(:,o) = mean(NSE(:,:,o),2);
end
for i=1:nrofi
    for o=1:nrofo
        FRFn(:,(i-1)*nrofo+o) = Yn(:,o)./Xs(:,i);   % noise frf calc
    end
end
                                  
end