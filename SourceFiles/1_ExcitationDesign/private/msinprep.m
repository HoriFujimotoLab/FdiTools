function [xtim,df]=msinprep(varargin)
%MSINPREP - Generate multisine time series from complex amplitudes.
%
% xtim      : generated time series
% df        : calculated common divider of frequencies
% Fdat      : output of msinclip containing freq points & Fourier amp
% N         : length of the time series (fs/df)
% fs        : sampling frequency (N*df)
% dev       : device for which the series is prepared:
%           'screen': plotting etc, no modification
%           'DAC'   : D/A converter, predistorted amp for inverse tf zoh
% 
% Usage: [xtim,df]=msinprep(Fdat,N,fs,dev);
% Example: [cx,crestx]=msinclip(fiddata([],[],[1:15]'/256));
%          arbitgen=msinprep(cx,512,1);
% Copyright (c) I. Kollar and VUBrussel, ELEC, 1991-2003
%
% See also: MSINCLIP, OPTEXCIT.


Fdat=varargin{1};
if isa(Fdat,'fiddata')|isstruct(Fdat) %new call
  error(nargchk(1,4,nargin))
  ch=Fdat.InputCharacter;
  if any(strmatch(ch,{'BL','Samples','Discrete'}))
    dev='screen'; %no predistrotion
  elseif strcmp(ch,'ZOH')|strcmp(ch,'FOH')
    dev='DAC';
  else
    dev='';
  end
  if (nargin>=4)&~isempty(varargin{4}), dev=varargin{4}; end
  if isempty(dev), dev='screen'; end
  if nargin<3, fs=[]; else fs=varargin{3}; end
  if nargin<2, N=[]; else N=varargin{2}; end
  cx=Fdat.Input;
  freqv=Fdat.FreqPoints;
else
  error(nargchk(2,5,nargin))
  if nargin<5, dev=''; else dev=varargin{5}; end
  if nargin<4, fs=[]; else fs=varargin{4}; end
  if nargin<3, N=[]; else N=varargin{3}; end
  cx=varargin{2};
  freqv=varargin{1};
  if isempty(dev), dev='DAC'; end
end
%
if length(freqv)~=length(cx), error('freqv and cx have different lengths'), end
[freqv,ind]=sort(freqv); freqv=freqv(:);
cx=cx(ind); cx=cx(:);
if any(freqv<0), error('freqv has negative element(s)'), end
freqvnz=freqv;
if freqvnz(1)==0, freqvnz(1)=[]; end
dfv=diff([0;freqvnz]);
if any(dfv==0), error('freqv must not contain equal elements'), end
%find common divider df
if isempty(N), remfnegl=max(freqv)/min(diff([0;freqvnz]))*eps*max(freqv);
elseif ~isempty(fs), remfnegl=fs/N/2;
else remfnegl=2*max(freqv)/N;
end
df=freqvnz(1);
if any(dfv<=0), error('freqv is not monotonously increasing'), end
while any(dfv>remfnegl)
  df0=df; df=min([df0;dfv]);
  dfv=sort(rem([df0;dfv],df));
  if isempty(N), remfnegl=max(freqv)/df*eps*max(freqv); end
  ind=find(dfv<=remfnegl); dfv(ind)=[];
end
fi=round(freqv/df); %harmonic numbers
if max(fi)>1023,
  warning(sprintf(['Maximum harmonic index found in ''msinprep''\n',...
            '   is %.0f, with df = %.3g Hz, T = %.3g s'],max(fi),df,1/df))
  if max(fi)>1e4
    disp('Large indexes are often due to inaccurately given frequency values.')
    disp('If this is the case, before invoking msinprep do')
    disp('   freqv = round(freqv*T)/T;')
    disp('where T is the desired period length.')
    disp(' ')
  end
end
devii=fi*df-freqv; ddf=fi\devii;
if max(abs(fi*df-freqv))>max(abs(fi*(df+ddf)-freqv)), df=df+ddf; end
maxdev=max(abs(fi*df-freqv));
if maxdev>2*eps*max(freqv)*max(fi)
  warning(sprintf('Maximum deviation of i*df from fi is %.2e in msinprep',maxdev))
end
%
if isempty(N)&isempty(fs) %none of them is given
  fs=df*2*(max(freqv)/df+1); N=round(fs/df);
end
if isempty(fs)
  fs=N*df;
elseif max(freqv)>(1+N*eps)*fs/2
  error('max of freqv exceeds fs/2')
end
if isempty(N)
  N=round(fs/df);
elseif N<(1-N*eps)*round(fs/df)
  error(sprintf('N (%.0f) must be at least %.0f for one period',...
      N,round(fs/df)))
end
N=round(N);
if abs(rem(df/(fs/N)+0.5,1)-0.5)>1e3*eps*N
  pno=N/(fs/df); digno=ceil(-log10(abs(rem(pno+0.5,1)-0.5)))+1;
  warning(sprintf(['The N (%.0f) samples cover %.',int2str(digno),'f'...
       ' periods,\n   this is not an integer number'],N,N/(fs/df)))
end
%
if strcmp(dev,'DAC') %precompensation for zero order hold
  tf=exp(-j*pi*freqv/fs).*sin(pi*freqv/fs+eps)./(pi*freqv/fs+eps);
  cx=cx./tf;
elseif strcmp(dev,'screen')|isempty(dev)
  %nothing to do
else
  error(['dev ''',dev,''' unknown'])
end
%
hno=round(freqv/(fs/N)); %harmonic numbers
Xk=zeros(N,1); Xk(hno+1)=2*cx; Xk(1)=Xk(1)/2; xtim=real(ifft(Xk));
if isa(Fdat,'fiddata')
  xtim=tiddata([],xtim,1/fs);
  if strcmp(dev,'screen'), xtim.inputcharacter='BL'; end
end
%%%%%%%%%%%%%%%%%%%%%%%% end of msinprep %%%%%%%%%%%%%%%%%%%%%%%%
