function [x,Xs,freqs,Xt,freqt] = msin(fs,df,fl,fh,itp,ctp,stp,Bn,An)
%MSIN - Multisine Excitation Signal generation.
% [x,Xs,freqs,Xt,freqt] = msin(fs,df,fl,fh,itp,ctp,stp,Bn,An)
% fs        : sample frequency in Hz
% df        : difference between spectral lines
% fl,fh     : minimal and maximal excitation frequency
% itp       : initial guess of phases (optional): 
%             's' is schroeder phases (default)
%             'r' is random phases between -pi and pi
% ctp       : 'c' = compressed multi sine (default)
%             'n' = not compressed multi sine
% stp       : type: 'o' = odd, 'f' = full multi sine (also default)
%             'O' = odd-odd, 'O2' = special odd-odd, 
% Bn,An     : continuous time transfer function model 
%             describing desired amplitude spectrum: | Bn / An |
% x,X       : time and frequency data of multisine
% freq      : excited frequency lines of multisine
% Xt,freqt  : total spectrum of multi sine
% Author    : Thomas Beauduin, KULeuven, 2015

% default values
if nargin<5, itp = 's'; end
if itp~='r', itp = 's'; end
if nargin<6, ctp = 'c'; end
if nargin<7, stp = 'full'; end
if nargin<8, Bn = 1; An = 1; end

% init calc
nrofs=fs/df;
n_max=round(fh/df)+1; 
n_min=ceil(fl/df)+1; 
w_s=1i*2*pi*df*(0:1:nrofs/2-1)';
w_s=w_s(1:n_max);
Mag = abs(polyval(Bn,w_s)./polyval(An,w_s));

% full multi sine flat spectrum
if stp == 'f'
   X0=zeros(n_max,1);
   X0(n_min:1:n_max) = Mag(n_min:1:n_max);
   ff = 1;
end;

% full multisine semilog spectrum (under-costruct)
if stp == 'l'
    X0=zeros(n_max,1);
end

% odd multi sine flat spectrum
if stp == 'o'
   vc = (1:2:n_max-1)';
   i = find(abs(vc-(n_min-1)) == min(abs(vc-(n_min-1))));
   i = max(i);
   ind = vc(i)+1;
   X0=zeros(n_max,1); 
   X0(ind:2:n_max) = Mag(ind:2:n_max); %ones(length(ind:2:n_max),1);
   ff = 2;   
end;

% odd-odd multi sine flat spectrum
if stp == 'O'
   vc = (1:4:n_max-1)';
   i = find(abs(vc-(n_min-1)) == min(abs(vc-(n_min-1))));
   i = max(i);
   ind = vc(i)+1;
   X0=zeros(n_max,1);
   X0(ind:4:n_max) = Mag(ind:4:n_max); %ones(length(ind:4:n_max),1);
   ff = 4;   
end;

% odd-odd multi sine flat spectrum
if stp == 'O2'
   vc1 = (1:8:n_max-1)';
   vc2 = (3:8:n_max-1)';
   vc = sort([vc1;vc2]);
   i = find(abs(vc-(n_min-1)) == min(abs(vc-(n_min-1))));
   i = max(i);
   ind = vc(i)+1;
   X0=zeros(n_max,1); %construction of the flat spectrum
   X0(ind:4:n_max) = Mag(ind:4:n_max);%ones(length(ind:4:n_max),1);
   ind = vc(i)+1;
   if max(vc(i)==vc1)      
      X0(ind:8:n_max) = Mag(ind:8:n_max);
      X0(ind+2:8:n_max) = Mag(ind+2:8:n_max);
   end;
   if max(vc(i)==vc2)
      X0(ind:8:n_max) = Mag(ind:8:n_max);
      X0(ind+6:8:n_max) = Mag(ind+6:8:n_max);
   end;
   ff = 4;
end;

% Phase Calculation
p=2; % starting value of p
X=X0;
if (ctp ~= 'n');
   while p<600          % Minimizes the l-2p norm
      X=msinl2pi(p,X,nrofs,[],0,[],10,1e-4,itp);
      p=ceil(p*2);      % Increment the value of p
   end
else
   if (itp=='r')
      X = randph(X);    % random phase 
   else
      X=schroed(X);     % If X=real then schroefrt
   end;
end;

% Total and Signal data
x = f2t(X,nrofs);
x=x/std(x);
Xt=t2f(x,nrofs);
Xs = Xt(n_min:n_max);
freqt=fs*(0:1:nrofs/2-1)'/nrofs; 
freqs=freqt(n_min:n_max);
