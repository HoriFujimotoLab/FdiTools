function [Yl,freql,Yo,freqo,Ye,freqe,Yn,freqn] = time2nld(x,y,fs,fl,fh,df)
%TIME2NLD - non-linear FRF distortions detection.
%
% x, y      : time data of odd-odd multisine measurement 
% fs,df     : sampling frequency and frequency resolution
% fl,fh     : lowest & highest frequency of excitated band
% Yl,freql  : linear system contribution to spectrum
% Yo,freqo  : odd non-linear contribution to spectrum
% Ye,freqe  : even non-linear contribution to spectrum
% Yn,freqn  : underlying noise contribution to spectrum
%

x=x(:); y=y(:);                                 % input/output vectoring
nrofs = fs/df;                                  % samples per period
nl = ceil(fl/df); nh = floor(fh/df);            % low & high freq points
freq = double((nl:1:nh)'/(nrofs/fs));           % freq lines
nr=length(freq);                                % number of freq lines
nrofp = double(floor(length(x)/nrofs));         % number of periods

% FFT Calculation
INPs=[]; OUTs=[]; OUT=[];
for i=1:nrofp
	INPs=[INPs fft(x(1+(i-1)*nrofs:i*nrofs))];  % fft of 1x period
	OUTs=[OUTs fft(y(1+(i-1)*nrofs:i*nrofs))];
end
for i=1:nrofp/2
    OUT=[OUT fft(y(1+(i-1)*nrofs*2:i*nrofs*2))];% fft of 2x period
end
INPs=INPs(nl+1:1:nh+1,:);                       % fft dc-term removal
OUTs=OUTs(nl+1:1:nh+1,:);
OUT=OUT(2*nl:1:2*nh+1,:);

m=1;
for k=1:nr*2
    if mod(k,2)~=0                              % uneven freq lines
        OUTn(m,:)=OUT(k,:);                     % noise data
        m=m+1;
    end
end
Xs=mean(INPs,2);                                % frequency averaging
Ys=mean(OUTs,2);
Yn=mean(OUTn,2);

% Contribution extraction
[~,s]=max(Xs(1:4));
switch s                                        % set initial freq 
    case 1, e = 1; o = 2;                           
    case 2, e =-1; o = 2;
    case 3, e =-1; o =-2;
    case 4, e =-3; o =-2;
end

freql=[]; Yl=[];                                % linear contribution
for i = s:4:length(freq) 
    freql = vertcat(freql,freq(i));
    Yl = vertcat(Yl,Ys(i)); 
end
freqe=[]; Ye=[];                                % even contribution
for i = s+e:2:length(freq)
    freqe = vertcat(freqe,freq(i));
    Ye = vertcat(Ye,Ys(i)); 
end
freqo=[]; Yo=[];                                % odd contribution
for i = s+o:4:length(freq)
    freqo = vertcat(freqo,freq(i));
    Yo = vertcat(Yo,Ys(i));
end
freqn = freq;                                   % underlying noise

end