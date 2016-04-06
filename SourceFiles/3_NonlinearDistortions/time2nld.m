function [Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,fs,fl,fh,df)
%TIME2NLD - non-linear FRF distortions detection (MIMO).
%   [Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,fs,fl,fh,df)
% x, y      : time data of odd-odd multisine measurement 
% fs,df     : sampling frequency and frequency resolution
% fl,fh     : lowest & highest frequency of excitated band
% Yl,freql  : linear system contribution to spectrum
% Yo,freqo  : odd non-linear contribution to spectrum
% Ye,freqe  : even non-linear contribution to spectrum
% Yn,freqn  : underlying noise contribution to spectrum
% Author    : Thomas Beauduin, KULeuven, PMA, 2014
[~,nrofi] = size(x);                                % number of inputs
[~,nrofo] = size(y);                                % number of outputs
nrofs = fs/df;                                      % samples per period
nl = ceil(fl/df); nh = floor(fh/df);                % low & high freq
freq = double((nl:1:nh)'/(nrofs/fs));               % full freq lines
nroff = length(freq);                               % number of freq lines
nrofp = double(floor(length(x)/nrofs));             % number of period

% Calculation of signal fft data
Xs = zeros(nroff,nrofi); Ys = zeros(nroff,nrofo);
INP = zeros(nroff,nrofp,nrofi); OUT = zeros(nroff,nrofp,nrofo);
for i=1:nrofi
    for p=1:nrofp
        Ip = fft(x(1+(p-1)*nrofs:p*nrofs,i));       % fft of 1x period
        INP(:,p,i) = Ip(nl+1:nh+1);                 % fft dc-term removal
    end
    Xs(:,i) = mean(INP(:,:,i),2);
end
for o=1:nrofo
    for p=1:nrofp
        Op = fft(y(1+(p-1)*nrofs:p*nrofs,o));       % fft of 1x period
        OUT(:,p,o) = Op(nl+1:nh+1);                 % fft dc-term removal
    end
    Ys(:,o) = mean(OUT(:,:,o),2);
end

% Calculation of noise fft data
Yn = zeros(nroff,nrofo);
OUT = zeros(nroff*2,floor(nrofp/2),nrofo);
NSE = zeros(nroff,floor(nrofp/2),nrofo);
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

% Extraction of non-linear contributions
[~,s] = max(Xs(1:4,1));
switch s                                            % set initial freq 
    case 1, even = 1; odd = 2;                           
    case 2, even =-1; odd = 2;
    case 3, even =-1; odd =-2;
    case 4, even =-3; odd =-2;
end

freqn = freq;                                       % underlying noise
freql=zeros(nroff/4,1); Yl=zeros(nroff/4,nrofo);    % lin contribution
freqe=zeros(nroff/2,1); Ye=zeros(nroff/2,nrofo);    % even contribution
freqo=zeros(nroff/4,1); Yo=zeros(nroff/4,nrofo);    % odd contribution
for o=1:nrofo
    index = 1;
    for f=s:4:nroff 
        freql(index) = freq(f); Yl(index,o) = Ys(f,o);
        index = index + 1;
    end
    index = 1;
    for f=s+even:2:nroff
        freqe(index) = freq(f); Ye(index,o) = Ys(f,o);
        index = index + 1;
    end
    index = 1;
    for f=s+odd:4:nroff
        freqo(index) = freq(f); Yo(index,o) = Ys(f,o);
        index = index + 1;
    end
end

end