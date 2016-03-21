function [cx,crx,crxmax,cry,crymax]=msinclip(varargin)
%MSINCLIP - Multisine design with mimimum crest factor.
%
%       Algorithm: time domain <--> frequency domain swapping, clipping in time
%       domain. The output is an object (structure) containing the complex Fourier
%       amplitudes; use MSINPREP to generate a tiime domain signal from these.
%
%       Output arguments:
%       cx = fiddata object of the designed multisine. In the property 'input'
%           it contains the coefficients of the complex Fourier series:
%           the absolute values are equal to halves of the real amplitudes).
%           For the ZOH case (see below for crestmode), precompensation
%           with the inverse of the ZOH transfer function is applied, and the
%           calculated coefficients are those of the ZOH-generated excitation
%           signal.
%       crinfo = information about the crest factors. Fields of this structure:
%         crx = crest factor of the generated multisine, calculated with the
%           given oversampling factor.
%         crxmax = worst case crest factor of the multisine
%         cry = crest factor of the multisine at the output of the linear
%           system, calculated with the given oversampling factor.
%         crymax = worst case output crest factor of the multisine
%
%       Input arguments:
%       Fdat = fiddata object (see 'help fiddata' and 'help(fiddata)',  
%           with the basic quantities (in stand-alone applications  
%           Fdat may be a structure with these fields):
%         freqpoints = vector of frequencies where the nonzero amplitudes are
%           given: the elements must be integer multiples of a df value, minimum
%           number of sines: 2. The frequency vector must be strictly
%           monotonously increasing.
%         input = absolute values of the desired nonzero complex amplitudes
%           (coefficients of the complex Fourier series at the corresponding
%           frequencies: halves of the real coefficients).
%           If any element is complex, the phases of the input vector will be
%           used as starting values, otherwise the Schroeder multisine is used.
%           For elements given with value NaN, auxiliary sine amplitudes are
%           returned ("snowing"): these are not included into the
%           calculation of the effective value of the useful excitation
%           signal, only help to decrease the peak value.
%           Default for these amplitudes: ones(length(fv),1)
%           Initial values for the snowing case can be given by a special
%           field of runmod: 'snowinginput' is the NaN-containing vector, the
%           real values are neglected.
%         output = (optional) outputs: input amplitudes multiplied by the complex
%           transfer function values at the given frequencies.
%           If the output is given, input-output optimization will be performed.
%       runmod = structure of run modifiers. Usually not necessary.
%         For the details, execute 'msinclip runmod'.
%
%       Usage: [cx,crinfo]=msinclip(Fdat,runmod);
%       Examples:
%         [cx,crestx]=msinclip(fiddata([],[],[1:15]'/256));
%         %Schroeder starting phases:
%         cx=msinclip(fiddata([],ones(12,1),4:15),struct('initset','Schroeder'));
%         %Prepare time series, 2048 points, fs = 1024 Hz:
%         timeobject=msinprep(cx,2048,1024);
%         %ZOH design:
%         cx=msinclip(fiddata([],ones(21,1)/sqrt(42),[0.2:0.01:0.4]*1e3),...
%              struct('crestmode','ZOH','N',100,'itno',350));
%       Copyright (c) I. Kollar and VUBrussel, ELEC, 1991-2004
%
%       See also: DIBS, MSINPREP, CRESTMIN.

if nargin==0, error('No input argument'), end
if (nargin==1)&isstr(varargin{1})
  if strncmpi(varargin{1},'preload',4)
    return %preload only
  elseif strcmp(varargin{1},'oldhelp')
    oldhelp msinclip
    return
  elseif strcmp(varargin{1},'runmod')|strcmp(varargin{1},'runpar')
    msincrmt
    return
  elseif strcmp(varargin{1},'secret')
    fprintf(['  ''Secret'' way of setting the number of iterations:\n',...
        '    global msinclip_max_iteration_number\n    msinclip_max_iteration_number=<iterations>;\n'])
    return
  end
  error(['First argument ''',varargin{1},''' is not allowed'])
end
%
global msinclip_max_iteration_number
c=computer; MatlV=version; isPC=strcmp(c(1:2),'PC'); isPlot=0; plotdone=0;
if isPC %Matlab5 under Win95 still greys out
  calldrawnow=1; calliterctrl=1;
else %regular case, e.g. Matlab5 on any platform, still may grey out
  calldrawnow=1; calliterctrl=1;
end
%calldrawnow=0; calliterctrl=0; %Force to eliminate drawnows
%calldrawnow=1; calliterctrl=1; %Force to call drawnows
%
itctrllastchecked='Continue'; termcond=0;
if calldrawnow, drawnow, end
%
crdef=1; %crdef=1: cr=max(abs(multisine))/effval, crdef=2:multisinepp/2/effval
%mbint(crdef) %assertion (neglected by Matlab 5.3, but not by 5.2)
yes=1; no=0;
global halffiginmsinclip
if ~exist('halffiginmsinclip'), halffiginmsinclip=''; end
if ~isempty(halffiginmsinclip)
  halffig=halffiginmsinclip; halffiginmsinclip='';
else halffig='n';
end
%for halffig='y', the plot of input minimization is subplot(1,2,1)
cldef='JS'; %clipping level; JS: Johan Schoukens, PG: Patrick Guillaume,
%                            VDO: Edwin van der Ouderaa
%
Fdat=varargin{1};
OK=1;
if isstruct(Fdat)
  if ~isfield(Fdat,'freqpoints'), OK=0; fprintf(['Fdat field ''freqpoints'' is missing'])
  elseif ~isfield(Fdat,'input'), OK=0; fprintf(['Fdat field ''input'' is missing'])
  elseif ~isfield(Fdat,'output'), Fdat.output='';
  end
end
if OK==0, error('First input argument is an improper structure'), end
initset='';
showmsgs=0;
Nmax=inf;
if exist('fdguidev.mat'), showmsgs=1; end %messages for development
if isa(Fdat,'fiddata')|isstruct(Fdat)
  if nargout>2, error('Too many output arguments'), end
  error(nargchk(1,2,nargin))
  initset='random'; %default for the new call
  fv=Fdat.freqpoints;
  ampv=Fdat.input;
  tf=Fdat.output;
  if ~isempty(tf)
    ind=find(~isnan(Fdat.input));
    tf(ind)=tf(ind)./Fdat.input(ind);
  end
  gmod=''; itno=[]; ovs=[]; N=[]; cl0=[]; sr=[]; pld=[]; crestmode='BL';
  if (nargin>=2)&~isempty(varargin{2})&isstruct(varargin{2})
    runmod=varargin{2};
    fn=fieldnames(runmod);
    for ii=1:length(fn)
      if strncmp(fn{ii},'initset',4)
        initset=getfield(runmod,fn{ii});
        if strncmpi(initset,'Schroeder',3)
        elseif strncmpi(initset,'random',4)
        elseif strncmpi(initset,'input',2)
        elseif isempty(initset)
        else error(['initset ''',initset,''' is invalid'])
        end
      elseif strncmp(fn{ii},'graphs',5), gmod=getfield(runmod,fn{ii});
      elseif strncmp(fn{ii},'showmsgs',5)
        showmsgs=getfield(runmod,fn{ii});
        showmsgs=~strncmp(showmsgs,'n',1);
      elseif strncmp(fn{ii},'crestmode',6)
        crestmode=getfield(runmod,fn{ii});
        if strcmpi(crestmode,'discrete'), gmod=[gmod,'ZOHd'];
        elseif strcmpi(crestmode,'ZOH')|strcmpi(crestmode,'DAC'), gmod=[gmod,'ZOHc'];
        elseif strncmpi(crestmode,'filtered-DAC',8)|strcmpi(crestmode,'fDAC')
          gmod=[gmod,'fDAC'];
        elseif strcmpi(crestmode,'BL') %default
        else error(['crestmode is not allowed: ''',crestmode,''''])
        end
      elseif strncmp(fn{ii},'itno',2), itno=getfield(runmod,fn{ii});
      elseif strncmp(fn{ii},'pmin',4) %for crestmin only
        warning('pmin has no meaning in msinclip')
      elseif strncmp(fn{ii},'relvarmax',4) %for crestmin only
        warning('relvarmax has no meaning in msinclip')
      elseif strncmp(fn{ii},'pmax',4) %for crestmin only
        warning('pmax has no meaning in msinclip')
      elseif strncmp(fn{ii},'oversampling',5), ovs=getfield(runmod,fn{ii});
      elseif strcmpi(fn{ii},'N'), N=getfield(runmod,fn{ii});
      elseif strcmpi(fn{ii},'Nmax'), Nmax=getfield(runmod,fn{ii});
      elseif strncmp(fn{ii},'cliplevel',5), cl0=getfield(runmod,fn{ii});
      elseif strncmp(fn{ii},'slewratemax',5), sr=getfield(runmod,fn{ii});
      elseif strncmp(fn{ii},'plottime',5), pld=getfield(runmod,fn{ii});
      elseif strncmp(fn{ii},'snowinginput',4)
        snowinp=getfield(runmod,fn{ii});
        if ~isempty(snowinp)
          if any(size(ampv)~=size(snowinp))
            error('snowinginput is inconsistent with the amplitudes')
          end
        end
        ind=find(isnan(snowinp))
        if ~isempty(snowinp), fv(ind)=-fv(ind); end
      else
        error(['Unrecognized runmod field ''',fn{ii},''''])
      end
    end
  end %for ii
  cry=[]; crymax=[];
else
  error(nargchk(1,10,nargin))
  if nargin<10, pld=[]; else pld=varargin{10}; end
  if nargin<9, sr=[]; else sr=varargin{9}; end
  if nargin<8, cl0=[]; else cl0=varargin{8}; end
  if nargin<7, N=[]; else N=varargin{7}; end
  if nargin<6, ovs=[]; else ovs=varargin{6}; end
  if nargin<5, itno=[]; else itno=varargin{5}; end %number of iterations
  if nargin<4, gmod=''; else gmod=varargin{4}; end
  if nargin<3, tf=[]; else tf=varargin{3}; end
  if nargin<2, ampv=[]; else ampv=varargin{2}; end
  fv=varargin{1};
end
if isempty(pld), pld=0; end
if isempty(sr), sr=inf; end
if sr<=0, error('sr is not positive'), end
if isnan(sr), sr=inf; end
if ~isempty(cl0)&(abs(cl0-.5)>=.5), error('Illegal cl0 value'), end
if isempty(N), N=2; end
if rem(N,2)~=0, error('N is not an even number'), end
if isempty(itno)|isnan(itno)
  if (length(msinclip_max_iteration_number)==1)&(msinclip_max_iteration_number>0)
    itno=msinclip_max_iteration_number;
  else
    itno=200;
  end
end
%
if isempty(gmod)|strcmp(gmod,'ZOHd')|strcmp(gmod,'ZOHc')|strcmp(gmod,'ZOH')|...
    strcmp(gmod,'BL')|strcmp(gmod,'fDAC')
  gmod=['graph10',gmod];
end
dtype='BL';
iZ=find(gmod=='Z');
if ~isempty(iZ)&(length(gmod)>=iZ+3)
  if strcmp(gmod(iZ+[0:3]),'ZOHd')
    dtype='ZOHd'; gmod(iZ+[0:3])=[]; iZ=[];
  end
end
if ~isempty(iZ)&(length(gmod)>=iZ+3)
  if strcmp(gmod(iZ+[0:3]),'ZOHc')
    dtype='ZOHc'; gmod(iZ+[0:3])=[]; iZ=[];
  end
end
if ~isempty(iZ)&(length(gmod)>=iZ+2)
  if strcmp(gmod(iZ+[0:2]),'ZOH')
    disp('In msinclip, ZOHd or ZOHc should be used instead of ZOH in gmod')
    %error(['gmod = ''',gmod,''' is not a valid option'])
    dtype='ZOHd'; gmod(iZ+[0:2])=[];
  end
end
if length(gmod)>=4
  indf=findstr(gmod,'fDAC');
  if ~isempty(indf), dtype='fDAC'; gmod(indf+[0:3])=[]; end
end
iB=find(gmod=='B');
if ~isempty(iB)&(length(gmod)>=iB+1)
  gmod(iB+[0:1])=[];
end
%
if strcmp(dtype,'BL')|strcmp(dtype,'fDAC')
  if isempty(ovs), ovs=16; end
else %ZOH
  if isempty(ovs), ovs=1; end
end
if ~strcmp(gmod,'nograph')&~strcmp(gmod,'graph')&~strcmp(gmod,'graph10')&...
        ~strcmp(gmod,'graph100')&~strcmp(gmod,'lastgraph')
  error(['gmod=''',gmod,''' is not a valid option'])
end
if isempty(ampv), ampv=ones(length(fv),1); end
%
if ovs<1, error('ovs must be at least 1'), end
if min(size(fv))>1, error('fv is not a vector'), end
indn=find(fv<0); %negative frequencies in signal start values for snowing
if ~isempty(indn)
  fv(indn)=-fv(indn);
  NaNv=NaN;
  cx=ampv; ampv(indn)=NaNv(ones(size(indn))); %simulate snowing input
  if strncmpi(initset,'Schroeder',3)|strncmpi(initset,'random',4), cx=[]; end %initset is stronger
else
  cx=[];
end
if min(size(ampv))>1, error('ampv is not a vector'), end
if length(fv)~=length(ampv)
  error('fv and ampv must have the same size')
end
fv=fv(:); ampv=ampv(:); tf=tf(:); %set to column vectors
if any(fv<0), error('fv has negative element(s)'), end
snow0='n'; freq0='n';
if fv(1)<10*realmin
  if ~isnan(ampv(1))
    dcx=real(ampv(1));
    if dcx==0, error('To define a zero dc value makes no sense'), end
  else
    dcx=0; dcy=0; snow0='y';
  end
  fv(1)=[]; ampv(1)=[]; if ~isempty(cx), cx(1)=[]; end
  if ~isempty(tf), dcy=dcx*tf(1); tf(1)=[]; end
  freq0='y';
else
  dcx=0; dcy=0;
end
%
ireg=find(~isnan(ampv));
isnow=find(isnan(ampv));
if isempty(ireg), error('No amplitude to optimize'), end
if sum(ampv(ireg)~=0)<1, error('No nonzero amplitude to optimize'), end
%
%find common divider df, not smaller than ovs*2*fmax/N if N is given.
df=dfcalc(fv);
fi=round(fv/df); %harmonic numbers
if max(fi)>1023,
  warning(sprintf(['Maximum harmonic index found in ''msinclip''\n',...
            '   is %.0f, with df = %.3g Hz, T = %.3g s'],max(fi),df,1/df))
  if max(fi)>1e4
    disp('Large indexes are often due to inaccurately given frequency values.')
    disp('If this is the case, before invoking msinclip do')
    disp('   fv = round(fv*T)/T;')
    disp('where T is the desired period length.')
    disp(' ')
  end
end
devii=fi*df-fv; ddf=fi\devii;
if max(abs(fi*df-fv))>max(abs(fi*(df+ddf)-fv)), df=df+ddf; end
maxdev=max(abs(fi*df-fv));
if maxdev>2*eps*max(fv)*max(fi)
  warning(sprintf('Maximum deviation of i*df from fi is %.2e in msinclip',maxdev))
end
if length(fi)<2, itno=0; end
if strcmp(dtype,'ZOHc')|strcmp(dtype,'ZOHd') %Zero-order hold
  if ovs==1, fsmod=1000*eps; else fsmod=0; end
  NZOH=2*ceil(ovs*max(fi)*(1+fsmod)); %FFT length
  Nov=NZOH;
else %Bandlimited multisine design or fDAC
  Nov=2^nextpow2(2*ovs*max([fi(ireg)+1;fi(isnow)])); %power of 2
end
Nt=max(N,Nov);
if strcmp(dtype,'BL')|strcmp(dtype,'fDAC')
  Nt=2^ceil(log2(Nt)-100*eps);
end
Nt=min(Nt,Nmax);
ovs=Nt/(2*max(fv)/df);
if isfinite(sr), srtxt=', with slew rate limitation'; else srtxt=''; end
if showmsgs==1
  if strcmp(dtype,'ZOHc')
    disp(['Continuous-time zero-order hold multisine design',srtxt])
  elseif strcmp(dtype,'ZOHd')
    disp(['Discrete-time zero-order hold multisine design',srtxt])
  elseif strcmp(dtype,'BL')
    disp(['Bandlimited multisine design',srtxt])
    Nt=2^ceil(log(Nt)/log(2)-100*eps);
  elseif strcmp(dtype,'fDAC')
    disp(['Multisine design for filtered DAC output',srtxt])
  end
  fprintf('df=%.5g, fmax=%.0f*df',df,max(fi(ireg)))
  if ~isempty(isnow), fprintf(', fmax(snow)=%.0f*df',max(fi(isnow))), end
  fprintf(', fsmin=%.0f*df\n',2*max([fi(ireg)+1;fi(isnow)]))
  fprintf('Prescribed oversampling: %.4g\n',ovs)
  if strcmp(dtype,'BL')|strcmp(dtype,'fDAC')
    fprintf(['Sample number rounded up to nearest power of 2,',...
        ' oversampling: %.3g\n'],ovs)
  end
  fprintf('FFT length in msinclip: %.0f\n',Nt)
end
dt=1/(Nt*df);
%
xeff=sqrt(2*sum(abs(ampv(ireg)).^2)+dcx^2);
srmin=sqrt(2*sum(abs(2*pi*fv(ireg).*ampv(ireg)).^2));
%Maximum overshoot between calculated values
if strcmp(dtype,'BL')|strcmp(dtype,'fDAC') %band-limited multisine design
  crdevx=sum(2*abs(ampv(ireg)).*(2*pi*fv(ireg)).^2*(dt/2)^2/2)/...
      sqrt(2*sum(abs(ampv(ireg)).^2));
  ge=' <= ';
else %zero-order hold
  crdevx=0; ge=': ';
end
if ~isempty(tf)
  yeff=sqrt(2*sum(abs(ampv(ireg).*tf(ireg)).^2)+dcy^2);
  crdevy=sum(2*abs(ampv(ireg).*tf(ireg)).*(2*pi*fv(ireg)).^2*(dt/2)^2/2)/...
    sqrt(2*sum(abs(ampv(ireg).*tf(ireg)).^2));
else
  crdevy=[];
end
%
%fiN is the pointer selecting the frequency points we are interested in
%fiv is a similar pointer for ampv
fiN=sort([fi;Nt-fi])+1;
fiNreg=sort([fi(ireg);Nt-fi(ireg)])+1;
fiNsnow=sort([fi(isnow);Nt-fi(isnow)])+1;
fiv=[[1:length(fi)]';[length(fi):-1:1]'];
fivreg=ireg([[1:length(fi(ireg))]';[length(fi(ireg)):-1:1]']);
fivsnow=isnow([[1:length(fi(isnow))]';[length(fi(isnow)):-1:1]']);
%
if isempty(cx) %initial values for snowing are not defined
  cx=zeros(size(ampv));
  if strncmpi(initset,'input',3) |...
      (isempty(initset) & (any(imag(ampv(ireg)))|any(real(ampv(ireg)))<0) )
    %starting phases are taken from input
    cx(ireg)=ampv(ireg);
  elseif strncmpi(initset,'random',4) %|isempty(initset)
    cx(ireg)=abs(ampv(ireg)).*exp(j*2*pi*rand(size(ampv(ireg))));
  else %definition of the Schroeder multisine
    pk=abs(ampv(ireg).^2)/sum(abs(ampv(ireg).^2));
    phi=zeros(length(fv(ireg)),1); %phi(1)=0
    for ni=2:length(fv(ireg))
      indl=find(ireg<ireg(ni));
      phi(ni)=phi(1)-...
        2*pi*sum(((fv(ireg(ni))-fv(ireg(1:ni-1)))/df).*pk(1:ni-1));
    end
    cx(ireg)=ampv(ireg).*exp(sqrt(-1)*phi);
  end
end
if strcmp(dtype,'ZOHc')|strcmp(dtype,'fDAC') %ZOH precompensation
  cx(ireg)=cx(ireg)./(exp(-j*pi*fv(ireg)*dt)...
        .*sin(pi*fv(ireg)*dt)./(pi*fv(ireg)*dt));
end
%
if calldrawnow, drawnow, end
if ~strcmp(gmod,'nograph')
  %First initialize plots
  if mean(get(gcf,'Color'))<=0.5, white='w'; else white='k'; end
  if get(0,'ScreenDepth')<4, blue=white; red=white; green=white;
  else blue='b'; red='r'; green='g';
  end
  %
  fsh=[fv';fv';fv']; fsh=fsh(:);
  ash=[zeros(1,length(fv));abs(ampv');zeros(1,length(fv))]; ash=ash(:);
  %Delete all axes objects in current figure
  %Use get because findobj may be missing (Matlab 4.1 or earlier)
  %delete(findobj(gcf,'Type','axes'));
  hax=get(gcf,'Children');
  if ~isempty(hax), for ih=1:length(hax)
    if strcmp(get(hax(ih),'Type'),'axes'), delete(hax(ih)), end
  end, end
  ithf=findobj('tag','msinclip_iterating'); delete(ithf)
  ithf=findobj('tag','fdident_dismiss'); delete(ithf)
  ithf=findobj('tag','signal_shape_dismiss_pb'); delete(ithf)
  if ~isempty(tf)
    subplot(2,1,1)
  end
end
incallow=12; %max. number of nonimproving steps before increase clipping level
nodecrease=0; %counter of nonimproving steps before increase clipping level
crxconst=0; %counter of no change steps
randtrial=0; %counter of random trials
randdone=0; %randomization just done
cly=0.99; %output clipping level set for the input only case
clxyuw1=0.1; clxyuw2=0.1; %clx and cly update weightings
if itno>0, stopiter=0; else stopiter=1; end %flag to terminate iteration
iopt=0; %number of optimal iteration cycle
im1=0; %-1 during stop by finish
lplt=clock; lplt(6)=lplt(6)-pld;
%initializations for the compiler:
pvx=[]; pvy=[]; txth=inf; crxold=0; cryold=0; dcrefx=0; dcrefy=0;
%
itcdone=0;
for i=[0:itno,itno] %the last cycle is for housekeeping and plotting the optimum
  %%
  if calliterctrl
    if strcmp(itctrllastchecked,'Finish')
      i=i-1-im1; im1=1; %correct last cycle
    end
  end
  if (i>=1) %for i=0 no clipping
    %set clipping level
    if strcmp(cldef,'VDO') %Method 1 for clipping level setting
      if nodecrease<incallow
        clx=interp1([0,1.6,2,3.5,4.5,1e6]',[.9,.9,.8,.3,.2,.2]',crx);
      else
        clx=min(clx+0.02,0.99);
      end
    elseif strcmp(cldef,'PG') %Method 2 for clipping level setting
      if i==1
        if ~isempty(cl0)
          clx=cl0;
        else
          clx=interp1([0,1.6,2,3.5,4.5,1e6]',[.9,.9,.8,.3,.2,.2]',crx);
        end
      end
      if crx<crxold %improvement found
        clx=clx/1.01;
      else
        clx=clx+clxyuw1*(1-clx);
      end
    elseif strcmp(cldef,'JS') %Method 3 for clipping level setting
      if i==1
        if ~isempty(cl0)
          clx=cl0;
        else
          clx=interp1([0,1.6,2,3.5,4.5,1e6]',[.9,.9,.8,.3,.2,.2]',crx);
        end
      end
      if crx<crxold %improvement found
        clx=max(0.2,clx-clxyuw2*(1-clx));
      else
        clx=clx+clxyuw1*(1-clx);
      end
    end
    if stopiter==0  %clip if still iterating
      if (~isempty(tf))|(randdone==1)
        %output step spoils multisine; otherwise it is already calculated
        randdone=0;
        %generate time function: multisine with useful effective value 1
        spect=zeros(Nt,1); spect(fiN)=cx(fiv);
        spect(Nt/2+1:Nt)=conj(spect(Nt/2+1:Nt));
        multisinex=real(ifft(spect))*Nt+dcx;
        %multisinex=multisinex/sqrt(mean(multisinex.^2));
        multisinex=multisinex/sqrt(2*sum(abs(cx(ireg)'*cx(ireg)))+dcx^2);
      end
      clms=multisinex; clind=find(abs(clms-dcrefx)>clx*pvx);
      minpts=5;
      if length(clind)>length(clms)-minpts,
        [absvalues,ind]=sort(abs(clms-dcrefx));
        clind=ind(minpts+1:length(ind)); clx=absvalues(minpts)/pvx;
      end
      clms(clind)=dcrefx+clx*pvx*sign(clms(clind)-dcrefx);
      if isfinite(sr) %slew rate limitation
        [srmaxi,imax]=max(xeff*abs(diff([clms;clms(1)]))/dt);
        srmax=srmaxi;
        Nsr=sum(abs(diff([clms;clms(1)]))/dt>sr/xeff);
        if i==1
          if showmsgs==1
            fprintf('Slew rate limitation, dt = %.4g, Amax = %.4g\n',...
              dt,xeff*max(abs(clms)))
            fprintf('  srdes = %.4g, lower bound on sr: %.4g\n',sr,srmin)
            if srmin>sr
              warning(sprintf('sr cannot be smaller than %.4g',srmin))
            end
          end
          fprintf(['  cyc: %.0f, srmax found: %.5g, number of points with ',...
              'large slew rate: %.0f\n'],i,srmax,Nsr)
        end
        Nldsp=0;
        while srmaxi>sr
          ind1=imax; ind2=ind1+1;
          if ind2>length(clms), ind2=1; end
          dclms=diff(clms([ind1,ind2]));
          clmod=0.999*0.5*sr/xeff*dt*sign(dclms);
          if Nsr<=10, clmod=clmod*0.9; end
          clms(ind1)=mean(clms([ind1,ind2]))-clmod;
          clms(ind2)=clms(ind1)+2*clmod;
          %Acceleration for many adjacent points
          if Nsr>10
            dind2l=[-Nsr+2:-1,1:Nsr-1]';
            ind2l=rem(length(clms)+ind1+dind2l-1,length(clms))+1;
            ind2lsr=find( sign(dclms)*abs(clms(ind2l)-clms(ind1))...
                              ./dind2l/dt > sr/xeff );
            if ~isempty(ind2lsr)
               clms(ind2l(ind2lsr))=clms(ind1)+...
                      0.999*sign(dclms)*sr/xeff*dt*dind2l(ind2lsr);
               %fprintf('  Nl = %.0f ',length(ind2lsr)); Nldsp=Nldsp+1;
            end
          end
          %
          [srmaxi,imax]=max(xeff*abs(diff([clms;clms(1)]))/dt);
        end %while
        if Nldsp>0, fprintf('\n'), end
        %fprintf('  srmaxi = %.4g\n',srmaxi)
      end %slew rate limitation
      %
      %find new complex amplitudes
      spect=fft(clms)/Nt;
      cx(isnow)=spect(fi(isnow)+1);
      %cx(ireg)=abs(cx(ireg))./abs(spect(fi(ireg)+1)).*spect(fi(ireg)+1);
      cx(ireg)=abs(cx(ireg)).*exp(j*angle(spect(fi(ireg)+1)));
      if strcmp(snow0,'y'), dcx=real(spect(1));
        if ~isempty(tf), dcy=dcx*tf(1); end
      end
    else
      if (i~=itno)&(termcond==1) %iteration number to be shown is the last one
        i=itno;
      end
    end %stopiter==0
  else %prepare for plot
    clx=1; crx=Nt;
  end
  %
  %generate time function: multisine with effective value 1
  spect=zeros(Nt,1); spect(fiN)=cx(fiv);
  spect(Nt/2+1:Nt)=conj(spect(Nt/2+1:Nt));
  multisinex=real(ifft(spect))*Nt+dcx;
  %multisinex=multisinex/sqrt(mean(multisinex.^2));
  multisinex=multisinex/sqrt(2*sum(abs(cx(ireg)'*cx(ireg)))+dcx^2);
  dcxn=mean(multisinex);
  %calculate crest factor
  minms=min(multisinex); maxms=max(multisinex);
  crxold=crx; %previous crest factor
  if crdef==1, %definition with reference to zero
    pvx=max(abs(minms),abs(maxms)); crx=pvx; dcrefx=0;
  else %crest factor is half of peak-to-peak value vs. effective value
    pvx=(maxms-minms)/2; crx=pvx; dcrefx=(minms+maxms)/2;
  end
  %administration of optimum
  if i==0 %still initialization
    cxopt=cx; crxopt=crx; dcxopt=dcx;
  elseif isempty(tf)
    if crx<crxopt-(1.01e-4)
      cxopt=cx; crxopt=crx; dcxopt=dcx; nodecrease=0; iopt=i;
    else
      nodecrease=nodecrease+1;
      if abs(crx-crxold)<10*eps, crxconst=crxconst+1; else crxconst=0; end
      %crxconst is the counter of successive crx steps with no change
    end
  end
  %
  %%
  if ~isempty(tf) %input-output iteration
    cy=tf.*cx;
    if i>=1 %for i=0 no clipping
      %set clipping level for output signal
      if strcmp(cldef,'VDO') %Method 1 for clipping level setting
        if nodecrease<incallow
          cly=interp1([0,1.6,2,3.5,4.5,1e6]',[.9,.9,.8,.3,.2,.2]',cry);
        else
          cly=min(cly+0.02,0.99);
        end
      elseif strcmp(cldef,'PG') %Method 2 for clipping level setting
        if i==1
          if ~isempty(cl0)
            cly=cl0;
          else
            cly=interp1([0,1.6,2,3.5,4.5,1e6]',[.9,.9,.8,.3,.2,.2]',cry);
          end
        end
        if cry<cryold
          cly=cly/1.01;
        else
          cly=cly+0.1*(1-cly);
        end
      elseif strcmp(cldef,'JS') %Method 3 for clipping level setting
        if i==1
          if ~isempty(cl0)
            cly=cl0;
          else
            cly=interp1([0,1.6,2,3.5,4.5,1e6]',[.9,.9,.8,.3,.2,.2]',cry);
          end
        end
        if cry<cryold %improvement found
          cly=max(0.2,cly-clxyuw2*(1-cly));
        else
          cly=cly+clxyuw1*(1-cly);
        end
      end %strcmp(cldef,...
      if stopiter==0  %clip if still iterating
        %generate time function: multisiney with effective value 1
        spect=zeros(Nt,1); spect(fiN)=cy(fiv);
        spect(Nt/2+1:Nt)=conj(spect(Nt/2+1:Nt));
        multisiney=real(ifft(spect))*Nt+dcy;
        %multisiney=multisiney/sqrt(mean(multisiney.^2));
        multisiney=multisiney/sqrt(2*sum(abs((cx(ireg).*tf(ireg))'...
                *(cx(ireg).*tf(ireg))))+dcy^2);
        clms=multisiney; clind=find(abs(clms-dcrefy)>cly*pvy);
        clms(clind)=dcrefy+cly*pvy*sign(clms(clind)-dcrefy);
        %find new complex amplitudes
        spect=fft(clms)/Nt;
        cy(isnow)=spect(fi(isnow)+1);
        cy(ireg)=abs(cy(ireg))./abs(spect(fi(ireg)+1)).*spect(fi(ireg)+1);
        cx=cy./tf;
        if strcmp(snow0,'y'), dcy=real(spect(1)); dcx=dcy/tf(1); end
      end
    else
      cly=1; cry=Nt;
    end
    %generate time function: multisiney with effective value 1
    spect=zeros(Nt,1); spect(fiN)=cy(fiv);
    spect(Nt/2+1:Nt)=conj(spect(Nt/2+1:Nt));
    multisiney=real(ifft(spect))*Nt+dcy;
    multisiney=multisiney/sqrt(mean(multisiney.^2));
    dcyn=mean(multisiney);
    %calculate crest factor
    minms=min(multisiney); maxms=max(multisiney);
    cryold=cry; %previous crest factor
    if crdef==1, %definition with reference to zero
      pvy=max(abs(minms),abs(maxms)); cry=pvy; dcrefy=0;
    else %crest factor is half of peak-to-peak value vs. effective value
      pvy=(maxms-minms)/2; cry=pvy; dcrefy=(minms+maxms)/2;
    end
    %administration of optimum
    if i==0 %still initialization
      cryopt=cry; nodecreasey=0;
    elseif max(crx,cry)<max(crxopt,cryopt)-(1.01e-4)
      %target: equal input and output crest factors
      cxopt=cx; dcxopt=dcx; crxopt=crx; cryopt=cry; nodecrease=0; iopt=i;
      if (strcmp(dtype,'BL')|strcmp(dtype,'fDAC'))&~isempty(isnow)
          crdevx=sum(2*abs(cx).*(2*pi*fv).^2*(dt/2)^2/2)/...
                        sqrt(2*sum(abs(cx).^2));
        if ~isempty(tf)
          crdevy=sum(2*abs(cx.*tf).*(2*pi*fv).^2*(dt/2)^2/2)/...
                sqrt(2*sum(abs(cx.*tf).^2));
        end
      end
    else
      nodecrease=nodecrease+1;
      if abs(crx-crxold)<10*eps, crxconst=crxconst+1; else crxconst=0; end
    end
  end %~isempty(tf)
  %
  %plot results in both cases
  if ((pld<=etime(clock,lplt)|(stopiter==1))&~strcmp(gmod,'nograph')) & ...
			( strcmp(gmod,'graph')|...
      (strcmp(gmod,'graph10')&((rem(i,10)==0)|(stopiter==1)))|...
      (strcmp(gmod,'graph100')&((rem(i,100)==0)|(stopiter==1)))|...
      (strcmp(gmod,'lastgraph')&(stopiter==1)) ) | isPlot
    lplt=clock;
    if strcmp(dtype,'ZOHc') %Zero-order hold
      tv=[0:Nt-1;1:Nt]/Nt/df; tv=tv(:);
      mv=[multisinex';multisinex']; mv=mv(:);
    else
      tv=[0:Nt]'/Nt/df; mv=[multisinex;multisinex(1)];
    end
    xlab=sprintf(['Eff. val.: %.4g, peak val.',ge,'%.4g'],...
            xeff,xeff*(crx+crdevx) );
    if isempty(tf)
      xlab=[xlab,sprintf(';  N=%.0f, ovsampl=%.2f',Nt,ovs)];
    end
    if ~isempty(tf)|strcmp(halffig,'y')|~isempty(isnow), subplot(2,1,1), end
    ha1=gca;
    hmscl=gcf;
    if itcdone==0 %try to put up iterctrl
      if ~any(findall(0,'type','uimenu','tag','iteration_menu'))&exist('iterctrl')
        iterctrl;
      end
      itcdone==1;
    end
    if ishandle(txth), delete(txth), end
    txth=axes('Position',[0.5,0,0.5,0.1],'parent',hmscl,'visible','off');
    ith=text(1,0,'Iterating...','parent',txth,'tag','msinclip_iterating',...
      'Verticalalignment','bottom','horizontalalignment', 'right');
    axes(ha1)
    %axes('Position',[0.13,0.5825,0.81,0.3175]); %upper subplot
    if ~strcmp(dtype,'ZOHd'), plot(tv,mv,['-',red])
    else plot(tv(1:length(tv)-1),mv(1:length(tv)-1),['+',red],...
              tv,mv,[':',red]) %discrete ZOH
    end
    %minms=min(mv); maxms=max(mv);
    grid off
    axv=axis; axv(2)=1/df;
    axv(4)=1.05*max(abs([minms,maxms])); axv(3)=-axv(4);
    axis(axv);
    xlabel(xlab)
    hold on
    plot([0,Nt-1]/Nt/df,dcrefx+clx*[pvx,pvx],[':',green],...
         [0,Nt-1]/Nt/df,dcrefx-clx*[pvx,pvx],[':',green],...
         [0,Nt-1]/Nt/df,[dcrefx,dcrefx],[':',white])
    hold off
    if strcmp(dtype,'BL')|strcmp(dtype,'fDAC') %BL multisine
      title([sprintf('Crf: %.5g (<=%.5g)',crx,crx+crdevx),...
           sprintf(', clip at: %.3g*max, iter: %.0f, opt. cyc: %.0f',...
              clx,i,iopt)])
    else %ZOH
      title([sprintf('Crf: %.5g',crx),...
           sprintf(', clip at: %.3g*max, iter: %.0f, opt. cyc: %.0f',...
              clx,i,iopt)])
    end
    ylabel('Normalized x(t)')
    if ~isempty(isnow)&isempty(tf) %plot snow amplitudes
      subplot(2,1,2)
      %plot([0;fsh;max(fsh)+df],[0;ash;0],['-',green])
      maxa=max(abs(cx));
      if ~exist('minlog'), minlog=[]; end
      if isempty(minlog), minlog=maxa/100; end
      if any(cx(isnow)~=0)
        minlog=min(minlog,max(abs(cx(isnow)))/10);
        minlog=max(minlog,maxa/10^5);
      end
      minlog=10^(floor(log10(minlog*1.001)));
      ind=find(ash<minlog); ash2=ash; ash2(ind)=minlog*ones(size(ind));
      semilogy([0;fsh;max(fsh)+df],[0;ash2;0],['-',green]), grid off
      axv=axis;
      axv(2)=max(fsh)+df;
      axv(3)=minlog;
      axv(4)=max(axv(4),(axv(4)/axv(3))^0.1*max(abs(ampv(ireg))));
      axis(axv);
      hold on
      if ~exist('fv2'), fv2=[]; end
      if isempty(fv2)
        fv2=[fv(isnow)';fv(isnow)';fv(isnow)';fv(isnow)']; fv2=fv2(:);
      end
      NaNv=NaN;
      asn2=[NaNv(:,ones(size(isnow')));minlog*ones(size(isnow'));...
              abs(cx(isnow))';NaNv(:,ones(size(isnow')))];
      asn2=asn2(:);
      semilogy(fv2,asn2,['-',red])
      hold off
         title(sprintf('Desired amplitudes: %.0f, snowing: %.0f',...
              length(ireg),length(isnow)))
      xlabel('Hz')
      if exist('fixsubplotlabels')
        fixsubplotlabels %make sure that title does not coincide with xlabel of other subplot
      end
    end
    if isempty(tf), drawnow, figure(gcf), pause(0), end
    plotdone=1;
  else
    plotdone=0;
  end %strcmp(gmod,'graph')
  if isempty(tf)&( strcmp(gmod,'lastgraph')|strcmp(gmod,'nograph')|...
              ( (stopiter==1)&(((clx==0.99)&(cly==0.99))|(i>=itno)) ) )
    %display value of crx
    if (rem(i,10)==0)|...
        ( (stopiter==1)...
        &(((clx==0.99)&(cly==0.99))|(i>=itno)) )
      if showmsgs==1
        fprintf('iteration=%.0f, crestx=%.5g',i,crxopt)
        if strcmp(dtype,'BL')|strcmp(dtype,'fDAC'), fprintf(', crestxmax=%.5g',crxopt+crdevx), end
        fprintf(', found in cycle %.0f, clipx=%.3g\n',iopt,clx)
      end
    end
  end
  %
  if ~isempty(tf)
    %plot results
		if ( ( (pld<=etime(clock,lplt))|(stopiter==1) ) & ~strcmp(gmod,'nograph') ) & ...
				(  strcmp(gmod,'graph')|...
				(strcmp(gmod,'graph10')&((rem(i,10)==0)|(stopiter==1)))|...
				(strcmp(gmod,'graph100')&((rem(i,100)==0)|(stopiter==1)))|...
				(strcmp(gmod,'lastgraph')&(stopiter==1))  ) ...
        | isPlot
      %
      xlab2=sprintf('Eff. val.: %.4g, peak val. <= %.4g',...
        yeff,yeff*(cry+crdevy) );
      xlab2=[xlab2,sprintf(';  Nt=%.0f, ovsampl=%.2f',Nt,ovs)];
      subplot(2,1,2), ha2=gca;
      if strcmp(dtype,'BL')|strcmp(dtype,'fDAC')
        plot([0:Nt]/Nt/df,[multisiney;multisiney(1)],['-',red])
      else
        plot([0:Nt-1]/Nt/df,multisiney,['+',red],...
                  [0:Nt]/Nt/df,[multisiney;multisiney(1)],[':',red])
      end
      %minms=min(multisiney); maxms=max(multisiney);
      grid off
      axv=axis; axv(2)=1/df;
      axv(4)=1.05*max(abs([minms,maxms])); axv(3)=-axv(4);
      axis(axv);
      xlabel(xlab2)
      hold on
      plot([0,Nt-1]/Nt/df,dcrefy+cly*[pvy,pvy],[':',green],...
           [0,Nt-1]/Nt/df,dcrefy-cly*[pvy,pvy],[':',green],...
           [0,Nt-1]/Nt/df,[dcrefy,dcrefy],[':',white])
      hold off
      title([sprintf('Output, crf: %.5g (<=%.5g)',cry,cry+crdevy),...
             sprintf(', clip at: %.3g*max',cly)])
      ylabel('Normalized y(t)')
      if exist('fixsubplotlabels')
        fixsubplotlabels %make sure that title does not coincide with xlabel of other subplot
      end
      drawnow, pause(0)
      zoom(gcf,'on')
    end
    if strcmp(gmod,'lastgraph')|strcmp(gmod,'nograph')|...
        ( (stopiter==1)&(((clx==0.99)&(cly==0.99))|(i>=itno)) )
      %display value of crx and cry
      if (rem(i,10)==0)|...
          ( (stopiter==1)&(((clx==0.99)&(cly==0.99))|(i>=itno)) )
        if showmsgs==1
          fprintf('iter=%.0f, crestx=%.5g, cresty=%.5g',i,crxopt,cryopt)
          fprintf(', found in cycle %.0f,\n   clipx=%.3g,',iopt,clx)
          fprintf(' clipy=%.3g',cly)
          if strcmp(dtype,'BL')|strcmp(dtype,'fDAC')
            fprintf(', crestxmax=%.5g, crestymax=%.5g\n',...
              crxopt+crdevx,cryopt+crdevy)
          else fprintf('\n')
          end
        end
      end
    end
  end %input-output optimization
  %
  if calliterctrl&~isempty(get(0,'Children'))
    if calldrawnow==1, drawnow, end
    %Store value to minimize iterctrl calls, to bypass bug on PC:
    if exist('iterctrl'), itctrlchecked=iterctrl('checked');
    else itctrlchecked='';
    end
    if strcmp(itctrlchecked,'Pause')|...
          strcmp(itctrlchecked,'Matlab prompt')
      if strcmp(itctrlchecked,'Pause'),
        disp('Iteration paused by GUI')
      else
        fprintf('Type your commands in the mini-window\n\n')
      end
      while strcmp(itctrlchecked,'Pause')|...
            strcmp(itctrlchecked,'Matlab prompt')
        drawnow %wait and draw
        if isPC
          %Such commands that allow to accept Enter in mini-window,
          %to bypass bug on PC
          clv=clock; drawnow, pause(0)
          while etime(clock,clv)<0.1, drawnow, pause(0), end
        end
        itctrlchecked=iterctrl('checked'); drawnow
      end %while
      disp('Continue...')
    end
    if strcmp(itctrlchecked,'Keyboard')
      disp('You may now enter commands within the workspace of msinclip')
      disp('Modify GUI selection and type ''return'' to continue')
      keyboard
      itctrlchecked=iterctrl('checked');
    end
    if strcmp(itctrlchecked,'Finish')
      %stopiter=1;
      randtrial=1; termcond=max(termcond,2);
      clx=0.99; cly=0.99; if i==0, itno=0; else itno=i; end
      itctrllastchecked='Finish'; iterctrl('Continue');
    elseif strcmp(itctrlchecked,'Cancel')
      iterctrl('Continue');
      disp('msinclip run canceled through GUI')
      cx=[]; crx=[]; crxmax=[]; cry=[]; crymax=[];
      if isa(Fdat,'fiddata')
        crinfo.crx=crx; crinfo.cry=cry; 
        crinfo.crxmax=crx+crdevx; crinfo.crymax=cry+crdevy;
        crx=crinfo;
      end
      return
    elseif strcmp(itctrlchecked,'Abort')
      iterctrl('Continue'); error('Abort requested through GUI')
    end
    isPlot=strcmp(itctrlchecked,'Plot');
  else
    isPlot=0;
  end
  %
  if calliterctrl&~isempty(get(0,'Children'))
    if calldrawnow==1, drawnow, end
    if strcmp(itctrlchecked,'Hold graph')&(plotdone==1)
      disp('Last graph of msinclip held by GUI')
      while strcmp(itctrlchecked,'Hold graph')
        drawnow
        itctrlchecked=iterctrl('checked');
      end %Wait
      disp('Continue...')
    end
    if strcmp(itctrlchecked,'Keyboard')
      disp('You may now enter commands within the workspace of elis')
      disp('Modify GUI selection and type ''return'' to continue')
      keyboard
      itctrlchecked=iterctrl('checked');
    end
    if strcmp(itctrlchecked,'Abort')
      iterctrl('Continue'); error('Abort requested through GUI')
    end
  end
  %
  if ((i~=0)&(clx>=0.99)&(cly>=0.99)) | (i>=itno) | (crxconst>10)
    %It is possible to stop iteration
    if (stopiter==1)&((i<itno)|(randtrial>0))
      %iteration converged, stop
      if ~strcmp(gmod,'nograph')
        xp=0.01; yp=0.01;
        graphh=gca;
        axes('Position',[0,0,1,1]); axis('off')
        if termcond==1
          text(xp,yp,'Iteration converged','VerticalAlignment','bottom')
        elseif termcond==2
          text(xp,yp,'Finish requested through GUI','VerticalAlignment','bottom')
        end
        axes(graphh)
      end
      if termcond==2, fprintf('Finish of msinclip requested through GUI\n'), end
      break
    else %reached termination criterion
      if (i<itno)&(randtrial<=0) %try random steps
        randtrial=randtrial+1; randdone=1;
        fprintf('%.0f. randomization done in msinclip in cycle %.0f\n',...
                         randtrial,i)
        stopiter=0; nodecrease=0; crxconst=0;
        cx=cxopt.*exp(sqrt(-1)*0.1*randn(length(cx),1));
        if strcmp(snow0,'y'), dcx=xeff/length(ampv)/10*2*(rand(1,1)-0.5); end
      else %prepare plot of best trial
        stopiter=1; cx=cxopt; dcx=dcxopt; crx=crxopt; crxmax=crx+crdevx;
        if ~isempty(tf), cry=cryopt; crymax=cry+crdevy; end
      end
    end
  end
  if itno==0
    if termcond==2, fprintf('Finish of msinclip requested through GUI\n'), end
    break
  end %exit if only plot was requested
end %i=1:itno
if strcmp(dtype,'ZOHc')|strcmp(dtype,'fDAC') %remove ZOH precompensation
  cx=cx.*(exp(-j*pi*fv*dt).*sin(pi*fv*dt)./(pi*fv*dt));
end
if strcmp(freq0,'y'), cx=[dcx;cx]; fv=[0;fv]; end
if calliterctrl
   hit=findobj(0,'tag','iteration_menu');
   if ~isempty(hit)&exist('iterctrl.m'), iterctrl('Continue'), end
end
%
if exist('ith') %Text 'Iterating' is shown
  if any(ishandle(ith))
    hp=get(ith,'parent'); delete(ith)
    set(hp,'units','pixels');
    pos=get(hp,'Position');
    p(4)=20; p(3)=50; p(2)=pos(2); p(1)=pos(1)+pos(3)-p(3);
    uicontrol('position',p,'string','Close','tag','fdident_dismiss',...
    'callback',...
    'delete(get(findobj(''tag'',''fdident_dismiss''),''parent'')'')')
  end
end
%
if isa(Fdat,'fiddata')|isstruct(Fdat) %new call
  if ~isempty(tf), y=cx.*tf; else y=[]; end
  if isstruct(Fdat),
    cx.Input=cx; cx.Output=y; cx.FreqPoints=fv;
  else
    cx=fiddata(y,cx,fv);
  end
  if strcmp(crestmode,'ZOH'), ich='ZOH';
  elseif strcmp(crestmode,'discrete'), ich='Discrete';
  else ich='BL';
  end
  cx.InputCharacter=ich;
  if ~isempty(tf)
    if strcmp(ich,'AASamples')|strcmp(ich,'BL'), och='BL';
    elseif strcmp(ich,'ZOH'), och='Samples';
    elseif strcmp(ich,'Discrete'), och='Discrete';
    end
    cx.OutputCharacter=och;
  end
  crinfo.crx=crx; crinfo.cry=cry;
  crinfo.crxmax=crx+crdevx; crinfo.crymax=cry+crdevy;
  crx=crinfo;
end
%%%%%%%%%%%%%%%%%%%%%%%% end of msinclip %%%%%%%%%%%%%%%%%%%%%%%%

function msincrmt
txt=str2mat('',...
'Possible fields and their values of the msinclip run modifier structure runmod',...
'==============================================================================',...
'initset = way of setting the inital values (overrides the real/complex',...
'  nature of the amplitudes).',...
'  Values: ''Schroeder'', ''random'' (default), ''input''',...
'graphs = (optional) if this is given with the value ''nograph'',',...
'  iteration results will not be plotted, if with the value',...
'  ''lastgraph'', the result of the last iteration only, if with',...
'  ''graph10'', of every 10th iteration, if with ''graph100'', every 100th ',...
'  iteration, if with ''graph'', a graph in every iteration will be',...
'  plotted. Default: ''graph10''',...
'showmsgs = show messages (''yes'' or ''no''), default: ''no''',...
'crestmode = way of crest factor minimization.',...
'  ''BL'': By default, a band-limited design is done.',...
'    We assume that the time domain signal will be generated from the',...
'    calculated amplitudes with very high clock frequency or with a',...
'    reconstruction filter assuming exact samples of a band-limited',...
'    signal. Thus, we assume that between the samples calculated in',...
'    msinclip the signal may be somewhat larger than the neighbouring',...
'    samples - this gives a limit for the crest factor - and no',...
'    predistortion is necessary.',...
'  ''ZOH'': we assume that a ZOH signal',...
'    will be generated from the samples, calculated with the same',...
'    clock frequency used in msinclip. In this case, no "overshoot"',...
'    can happen between the samples, but the amplitudes will be',...
'    distorted by the ZOH, so precompensation is introduced for this,',...
'    using the reciprocal of  tf_ZOH = sin(pi*fv*dt)./(pi*fv*dt).',...
'  ''filtered-DAC'' or ''fDAC'': we assume that a staircase signal',...
'    will be generated from the samples, calculated with the same',...
'    clock frequency used in msinclip, and then lowpass filtering',...
'    is applied. In this case, the signal is generated using a DAC,',...
'    but, because of the lowpass filter applied afterwards, "overshoot"',...
'    can happen between the samples. Thus, the amplitudes will be',...
'    distorted by the ZOH, so precompensation is introduced for this,',...
'    using the reciprocal of  tf_ZOH = sin(pi*fv*dt)./(pi*fv*dt).',...
'  ''discrete'': we deal with a discrete signal.',...
'    No predistortion is necessary, and no overshoot may happen.',...
'itno = (optional) number of iteration cycles, default: 200. If itno=0,',...
'  the starting values will be returned (Schroeder or external) ',...
'oversampling = (optional) minimal oversampling factor: the sampling frequency',...
'  will be   fs>=ovs*2*fmax,  default: 16 for BL design, 1 for ZOH',...
'N = number of points in the time series (even number); for band-limited',...
'  design this is the minimum number of points to perform base-2 FFT',...
'  effectively',...
'Nmax = maximum number of points in the time series (even number)',...
'cliplevel = (optional) initial clipping level, 0<cl0<1',...
'slewratemax = (optional) allowed maximum of the slew rate (maximum ',...
'  steepness of the waveform). After each clipping, the routine limits the',...
'  slew rate of the signal to the given maximum value. This does not ',...
'  guarantee that the slew rate stays indeed below the limit after ',...
'  restoration of the amplitudes, but usually drives the design towards ',...
'  a lower slew rate than it would be without slew rate limitation.',...
'plottime = minimum time between plots in seconds',...
'snowinginput = a vector which marks the snowing amplitudes as NaN''s');
disp(txt)

function [df,fi]=dfcalc(freqvect,digits,Nfs)
%DFCALC Calculates maximum common divisor in frequency vector
%
%       [df,fi]=DFCALC(freqvect,digits,Nfs)
%
%       Output arguments:
%       df = maximum common divisor (reciprocal of period length)
%       fi = frequency indices
%
%       Input arguments:
%       freqvect = vector of frequencies
%       digits = number of exact digits in the data
%       Nfs = maximum harmonic number associated to the sampling frequency
%             (which is larger than two times the maximum frequency)
%             The algorithm only makes an attempt, but does not guarantee to
%             fulfill the prescription
%
%       Usage: [df,fi]=dfcalc(freqvect,digits,Nfs);
%       Example:
%         df=dfcalc([5.1:2:16]);
%
%       See also: MSINCLIP.

%       Copyright (c) I. Kollar and Vrije Universiteit Brussel, ELEC, 1991-2001
%       All rights reserved.
%       $Revision: $
%       Last modified: 05-Jul-2001

if nargin<3, Nfs=[]; end, if isempty(Nfs), Nfs=inf; end
if Nfs<=2, error('Nfs<=2'), end
if nargin<2, digits=[]; end, if isempty(digits), digits=17; end
if iscell(freqvect)
  fsave=freqvect; freqvect=[];
  for ii=1:prod(size(fsave))
    freqvect=[freqvect;fsave{ii}];
  end
end
if isempty(freqvect), df=[]; return, end
if min(size(freqvect))~=1, error('freqvect is not a vector'), end
if length(freqvect)==1, df=freqvect; fi=1; return, end
freqvect=sort(freqvect);

%find common divider df, not smaller than 2*fmax/(Nfs-2) if Nfs is given.
%remfnegl is the remainder which is negligible in the frequency vector
emax=10^(-digits)*max(freqvect);
[mfv,mfvind]=min(freqvect);
if mfv==0, freqvect(mfvind)=[]; end
dfreqvect=diff(sort([0;freqvect(:)]));
ind=find(dfreqvect==0); if ~isempty(ind), dfreqvect(ind)=[]; end
df=min(dfreqvect);
if (Nfs==2)|~isfinite(Nfs),
  remfnegl=max(freqvect)/df*eps*max(freqvect)*length(freqvect);
else
  remfnegl=2*max(freqvect)/(Nfs-2)/2;
end
while any(dfreqvect>remfnegl)
  df0=df; df=min([df0;dfreqvect]);
  if (Nfs==2)|~isfinite(Nfs),
    remfnegl=max(freqvect)/df*eps*max(freqvect)*length(freqvect);
  end
  remfnegl=max(remfnegl,emax);
  dfreqvect=sort(abs(rem([df0;dfreqvect]+df/2,df)-df/2));
  ind=find(dfreqvect<=remfnegl); dfreqvect(ind)=[];
end
%
fi=round(freqvect/df); %harmonic numbers
ind=find(fi~=0);
df=mean(freqvect(ind)./fi(ind));
if max(fi)>=1e5
  disp(sprintf('Warning! the maximum harmonic number is %.0f.',max(fi)))
  disp('Large indexes are often due to inaccurately given frequency values.')
  disp('If this is the case, before invoking dfcalc round the frequencies:')
  disp('   freqvect = round(freqvect*T)/T;')
  disp('where T is the desired period length.')
  disp(' ')
end


function message = oldhelp(mfile,opt)
%OLDHELP MFILE  displays the help on the old call form of fdident functions

%       Copyright (c) I. Kollar and Vrije Universiteit Brussel, ELEC, 1996-2000
%       All rights reserved.
%       $Revision: $
%       Last modified: 09-Dec-2000

if nargin==0, mfile=''; end
if nargin<2, opt=''; end
if isempty(mfile)
  help oldhelp
  return
elseif isempty(opt)&any(strmatch(mfile,...
    {'crestmin'}))
  %Just <funcname> oldhelp
  feval(mfile,'oldhelp')
  return
end
if isstr(mfile)&(length(mfile)>0)&isstr(opt)&(length(opt)>0)&~strcmpi(opt,'date')
  %help with option
  eval([mfile,' ',opt])
  return
end
mfilein=mfile;
perpos=find(mfile=='.');
if isempty(perpos), mfile=[mfile,'.m']; end
perpos=find(mfile=='.');
mfname=mfile(1:perpos(1)-1);
clear perpos
%
mfline=which(mfile);
if isempty(mfline)&isdir(mfilein)
  if strcmp(mfilein(end),filesep)|strcmp(mfilein(end),'\')|...
      strcmp(mfilein(end),'/')
    %mfilein(end)='';
  else
    mfilein=[mfilein,filsesep];
  end
  type([mfilein,'Contents.m'])
  return
end
if isempty(findstr(['fdident',filesep,'fdident',filesep,lower(mfile)],lower(mfline))) &...
    isempty(findstr(['fdident',filesep,'fddemos',filesep,lower(mfile)],lower(mfline))) &...
    isempty(findstr([filesep,'fd',filesep,lower(mfile)],lower(mfline))) &...
    isempty(findstr([filesep,'fdd',filesep,lower(mfile)],lower(mfline))) & ...
    ~strcmpi(opt,'date')
  %if isempty(mfline), error(['Not fdident M-file: ',mfile]), end
  %error(['Not fdident M-file: ',mfline])
end
%
[fid,errmess]=fopen(mfile,'r');
if fid==-1
  if exist(mfname)==5
    messagei=[mfname,' is a built-in function, no M-file available'];
    disp(messagei)
    if nargout>0, message=messagei; end
    help(mfile)
    return
  elseif exist(mfname)==1
    if strcmp(mfname,'mfile')|strcmp(mfname,'mfname')|...
        strcmp(mfname,'nargin')|strcmp(mfname,'nargout')|...
        strcmp(mfname,'fid')|strcmp(mfname,'errmess')
      messagei=['File ',mfile,': ',errmess];
    else
      messagei=[mfname,' is a built-in variable, no help M-file available'];
    end
    disp(messagei)
    if nargout>0, message=messagei; end
    return
  elseif exist(mfname)==0
    c=computer;
    pathsep = ':'; %Unixpath separator character
    dirsep = '/'; %Unix directory separator character
    if strcmp(c(1:2),'PC')
      pathsep = ';'; dirsep = '\';
    elseif strcmp(c(1:2),'MA')
      pathsep = ';'; dirsep = ':';
    end
    ps=[pathsep,path,pathsep];
    ind=find(ps==pathsep);
    dirfound=0;
    for i=2:length(ind)
      if ind(i)>length(mfname)
        if strcmp(ps(ind(i)-[length(mfname):-1:1]),mfname)
          dirfound=1;
          break
        end
      end
    end %for i
    if dirfound==1
      messagei=[mfname,' is not an M-file, it is a directory:',sprintf('\n'),...
          ps(ind(i-1)+1:ind(i)-1)];
    else
      messagei=[mfname,' does not exist in workspace of function OLDHELP'];
    end
    disp(messagei)
    if nargout>0, message=messagei; end
    return
  end
  disp(errmess)
  if nargout>0, message=str2mat(message,errmess); end
  return
end
%
filestr=fread(fid); fclose(fid);
filestr=setstr(filestr');
lfstr=length(filestr);
if lfstr==0, error(['File ',mfile,' is empty']), end
if strcmpi(mfname,'elis'), scanlen=14000;
else scanlen=7000;
end
filestrsh=filestr(1:min(scanlen,length(filestr)));
indcr=[find((filestrsh==13)|(filestrsh==10)),length(filestrsh)+1]; %cr or lf
helpbegind=1; %index of first character of help to scan for call forms
helpendind=0; %index of last character of help to scan for call forms
if strcmpi(opt,'date')
  %Look for the file 'Last modified:'
  ind=findstr('last modified:',lower(filestrsh)); fmod=15;
  if length(ind)>1, ind=ind(1); end
  if isempty(ind), ind=findstr('modified:',lower(filestrsh)); fmod=10; end
  if length(ind)>1, ind=ind(1); end
  if isempty(ind), ind=findstr('date:',lower(filestrsh)); fmod=6; end
  if length(ind)>1, ind=ind(1); end
  if ~isempty(ind)
    %ind=ind+fmod;
    indind=min(find(indcr>ind));
    inde=indcr(indind)-1;
    dstr=filestrsh(ind:inde);
    indd=find(dstr=='$'); if ~isempty(indd), dstr(indd)=''; end
    disp(dstr);
  end
  return
end
%
%We will look for the word '%' and 'function' or 'Old fdident help' as string
%We will also explore here the beginning and end of the help text
firstchar=findstr(filestrsh,['%','function']); %Different old version definition
if isempty(firstchar), firstchar=findstr(filestrsh,['%','Old fdident help']); end
nfdstr=['%','Old fdident help: same as now'];
if ~isempty(firstchar)&strcmp(filestrsh(firstchar(1)+[0:length(nfdstr)-1]),nfdstr)
  help(mfile), return
elseif ~isempty(firstchar)
  firstchar=firstchar(1);
else
  %disp(' ')
  %disp('Warning: Modified fdident help not found, regular help is invoked')
  help(mfile), return
end
%
fhc=min(find(filestrsh(firstchar+2:length(filestrsh))=='%'));
nextchi=firstchar+1+fhc; %first 'hidden' help character: % sign
fhc=nextchi;
%
nch=filestrsh(nextchi); cri=1; stopsig=0; inhelp=1; funfound=-1;
commentline=1;
while stopsig==0
  if ( (nch==' ') | ((nch>=9)&(nch<=13)) ) %blank character
    commentline=commentline-1;
    if (inhelp==1)&(commentline<=0), inhelp=0; helpendind=nextchi-1; end
    nextchi=nextchi+1;
    if nextchi<=lfstr
      if (nch==13)&(filestrsh(nextchi)==10) %cr-lf on PC
        nextchi=nextchi+1;
      end
    end
  elseif nch=='%'
    if helpbegind==1, helpbegind=nextchi+1; inhelp=1; end
    nextchi=indcr(min(find(indcr>nextchi)));
    commentline=2;
  else %general character, cycle has to stop
    stopsig=1;
  end
  %Prepare next cycle
  if (nextchi>length(filestrsh)-400)&(length(filestrsh)<length(filestr))
    filestrsh=filestr;
    indcr=[find((filestrsh==13)|(filestrsh==10)),length(filestrsh)+1]; %cr or lf
  end
  if nextchi<=lfstr
    nch=filestrsh(nextchi);
  else
    stopsig=1;
  end
end %while
if inhelp==1, helpendind=nextchi-1; end
%
helptext=filestrsh(fhc+1:helpendind);
if length(helptext)>1;
  ind=findstr(helptext,setstr([13,10]));
  for ii=length(ind):-1:1, helptext(ind(ii))=''; end
end
ind=[findstr(helptext,[setstr(13),'%']),findstr(helptext,[setstr(10),'%'])];
for ii=length(ind):-1:1, helptext(ind(ii)+[0,1])='\n'; end
ind=findstr(helptext,'%');
for ii=length(ind):-1:1
  helptext=[helptext(1:ind(ii)),helptext(ind(ii):length(helptext))];
end
fprintf(['\n',helptext,'\n'])
%
ci=findstr('Copyright',filestrsh);
bci=max(find(filestrsh(1:ci(1))=='%')); %begin index of copyright string
yi=min(findstr('odified',filestrsh(bci+1:length(filestrsh))))+bci;
if isempty(bci), return, end
eci=min([find(filestrsh(yi+[0:min(length(filestrsh)-yi,80)])==setstr(13)),...
    find(filestrsh(yi+[0:min(length(filestrsh)-yi,80)])==setstr(10))])+yi-2;
ctext=filestrsh(bci:eci);
%
ctextm=ctext;
%Fix strange bug
eol=find((ctext==setstr(13))|(ctext==setstr(10)));
for ii=1:length(eol)
  if ii==1, ctextm=ctext(1:eol(1)-1); cti=eol(1)+1;
  else
    if eol(ii)>eol(ii-1)+1
      ctextm=[ctextm,sprintf('\n'),ctext(cti:eol(ii)-1)]; cti=eol(ii)+1;
    else
      cti=cti+1;
    end
  end
end %for
if cti<=length(ctext), ctextm=[ctextm,sprintf('\n'),ctext(cti:length(ctext))]; end
ind=find(ctextm=='%');
ctextm(ind)='';
%disp(ctextm) %display copyright info
%
%End of oldhelp

%%%%%%%%%%%%%%%%%%% End of file msinclip.m %%%%%%%%%%%%%%%%%%%%%%%%
