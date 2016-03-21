function X=msinl2pi(p,X,Nmax,Fa,W,H,iterno,relvar,itp)
%MSINL2PI - Lp-norm optimization of multisine phase.
%
% p         : Value of p in L2p
% X         : Matrix with input spectra (1=DC)
% Nmax      : Max. number of points in time domain
% Fa        : Set of the additional harmonic numbers (DC = 1)
% W         : Weighting row vector, e.g., W=[1=same_cr, 2=same_a];
% H         : Transfer function matrix 
% iterno    : Max. number of iterations
% relvar    : Max. relative variation of L2p-norm
% itp       : initial guess of phases (optional)
% X         : Matrices with optimal input spectra
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%
% See also SCHROED, RANPH, CREST.
%%%%%
if nargin<8, relvar=1e-6; end      % Default rel. dev. cost f.
if nargin<7, iterno=10; end        % 10 iterations by default
if nargin<6, H=[]; end             % Only input cf-minimization
if nargin<5, W=0; end, W=W(:)';    % W is a row vector (or scalar)
if nargin<4, Fa=[]; end, Fa=Fa(:); % Fa is a column vector
if nargin<3, Nmax=2048; end
if nargin<9, itp = 's'; end
if itp~='r', itp = 's'; end
%       Fe : Set of the effective harmonic numbers (index 1 = DC)
[rowno,colno]=size(X);  % Fe = all indices of non-zero X's entries - Fa
if colno==1, dummy=X; else dummy=max(abs(X)')'; end
Fe=cumsum(ones(size(dummy)));
Fe=Fe(dummy~=0);
clear dummy
for I=1:length(Fa), Fe(Fe==Fa(I))=[]; end  % Elimination of Fa indices
Feno=length(Fe); % Number of effective frequencies
Fano=length(Fa); % Number of additional frequencies
Ft=[Fe;Fa];      % Set of all the considered spectral lines
Ftno=Feno+Fano;  % Total number of considered spectral lines
kmax=max([max(Fe),max(Fa)])-1;  % Largest harmonic number of X
if length(X(:,1))<=kmax, X(kmax+1,1)=0; end
Ndummy=2*p*kmax+1; % Number of points needed for identity between 
     % L2p (summation) and l2p (integration) criterions
N=2; while Ndummy>=N, N=N*2; end % Select N=2^x > Ndummy
N=min([N,Nmax]);                 % If N>Nmax then N:=Nmax
if all(imag(X)==0), 
  if (itp=='r') 
    X = randph(X); % If X=real then randph
  else
    X=schroed(X); % If X=real then schroed
  end;
end;
   
if ~isempty(H)
   H=H(1:kmax+1);% Correct range setting of H
   Y=zeros(size(X));
   Y=H.*X;           % Output spectra
   if W==0, W=[effval(X,Fe),effval(Y,Fe)].^(-1); end
   % Same crest factors
   if W==1, W=[effval(X,Fe),effval(X,Fe)].^(-1); end
   % Same extreme value
   X=X*W(1);Y=Y*W(2);
   x=f2t(X,N);               % Input signals
   y=f2t(Y,N);               % Output signals
else
   W=1/effval(X,Fe);
   X=X*W;        % Scaling so that max|x(t)|~1
   x=f2t(X,N);               % Input signals
end

if isempty(H)
   cost=lpnorm(x(:),2*p);        % Calculation of the L2p-norm
else
   cost=lpnorm([x(:);y(:)],2*p);
end
                                % (L2p=costfunction)
X0=X;cost0=cost;iter0=0;        % Initialization of the best results
                                % (Best = lowest cost function)
indexmat=(Ft-1)*ones(size(Ft))';% Initialization of the (k-l)-matrix
indexmin=indexmat-indexmat'+N;  % and the (k+l)-matrix (Offset=N)
indexplus=indexmat+indexmat'+N;
phkphlmin=indexmin(:);          % Transformation of the matrices
phkphlplus=indexplus(:);        % into a vector
if ~isempty(Fa)
   phkalmin=indexmin(:,Feno+1:Ftno);
   phkalplus=indexplus(:,Feno+1:Ftno);
   akalmin=indexmin(Feno+1:Ftno,Feno+1:Ftno);
   akalplus=indexplus(Feno+1:Ftno,Feno+1:Ftno);
   phkalmin=phkalmin(:);        % Transformation of the matrices
   phkalplus=phkalplus(:);      % into a vector
   akalmin=akalmin(:);
   akalplus=akalplus(:);   
end
relax=0.01;skip=0;              % Levenberg-Marquardt parameters

iter=0;                         % Iteration counter

relerror=inf;                   % Relative deviation of the cost function

while (iter<iterno)&&(relerror>relvar)  % MAIN ITERATION LOOP
   iter=iter+1;
   X2p2=fft(x.^(2*p-2));               % FFT of x^(2p-2)
   X2p2=[flipud(conj(X2p2(2:N)));X2p2];
   dummyt=X(Ft);
   
   %  Calculation of JphJph and JphE
   Qmin=conj(dummyt*dummyt');
   Qplus=conj(dummyt*dummyt.');
   Qmin(:)=Qmin(:).*X2p2(phkphlmin);
   Qplus(:)=Qplus(:).*X2p2(phkphlplus);
   JphJph=p*real(Qmin-Qplus);
   JphE=sum(imag(Qmin+Qplus)')';
   if ~isempty(Fa)
      dummya=exp(1i*angle(X(Fa)));
      
      %  Calculation of JphJa and JaE
      Qmin=conj(dummyt*dummya');
      Qplus=conj(dummyt*dummya.');
      Qmin(:)=Qmin(:).*X2p2(phkalmin);
      Qplus(:)=Qplus(:).*X2p2(phkalplus);
      JphJa=p*imag(Qmin+Qplus);
      JaE=sum(real(Qmin+Qplus))';
      %  Calculation of JaJa
      Qmin=conj(dummya*dummya');
      Qplus=conj(dummya*dummya.');
   Qmin(:)=Qmin(:).*X2p2(akalmin);
      Qplus(:)=Qplus(:).*X2p2(akalplus);
      JaJa=p*real(Qmin+Qplus);
   end

   if ~isempty(H)
      Y=(H.*X)*W(2)/W(1);
      y=f2t(Y,N);
      Y2p2=fft(y.^(2*p-2));            % FFT of y^(2p-2)
      Y2p2=[flipud(conj(Y2p2(2:N)));Y2p2];
      dummyt=Y(Ft);

      %  Calculation of JphJph and JphE
      Qmin=conj(dummyt*dummyt');
      Qplus=conj(dummyt*dummyt.');
      Qmin(:)=Qmin(:).*Y2p2(phkphlmin)
      Qplus(:)=Qplus(:).*Y2p2(phkphlplus);
      JphJph=JphJph+p*real(Qmin-Qplus);
      JphE=JphE+sum(imag(Qmin+Qplus)')';
      
      if ~isempty(Fa)
        dummya=exp(1i*angle(Y(Fa)));
        
        %  Calculation of JphJa and JaE
        Qmin=conj(dummyt*dummya');
        Qplus=conj(dummyt*dummya.');
        Qmin(:)=Qmin(:).*Y2p2(phkalmin);
        Qplus(:)=Qplus(:).*Y2p2(phkalplus);
        JphJa=JphJa+p*imag(Qmin+Qplus);
        JaE=JaE+sum(real(Qmin+Qplus))';

        %  Calculation of JaJa
        Qmin=conj(dummya*dummya');
        Qplus=conj(dummya*dummya.');
        Qmin(:)=Qmin(:).*Y2p2(akalmin);
        Qplus(:)=Qplus(:).*Y2p2(akalplus);
        JaJa=JaJa+p*real(Qmin+Qplus);
      end
   end

   if isempty(Fa)
      A=JphJph;b=JphE;
   else
      A=[JphJph,JphJa;JphJa',JaJa];b=[JphE;JaE];
   end

   diagA=diag(A);       
   A=A+relax*diag(diagA+max(diagA)*eps); % Adding diagonal matrix to A

   % (Levenberg-Marquardt)
   Delta=-A\b;
   X(Fe)=abs(X(Fe)).*exp(1i*(angle(X(Fe))+Delta(1:Feno))); % Updating of X
   X(Fa)=(abs(X(Fa))+Delta(Ftno+1:Ftno+Fano))...
       .*exp(1i*(angle(X(Fa))+Delta(Feno+1:Ftno)));
   if ~isempty(H),Y=(H.*X)*W(2)/W(1);y=f2t(Y,N);end     % Updating of Y,y
   x=f2t(X,N);                                          % Updating of x
   
   if isempty(H)
      cost=lpnorm(x(:),2*p);        % Calculation of the L2p-norm
   else
      cost=lpnorm([x(:);y(:)],2*p);
   end
   
   relerror=abs(cost-cost0)/cost0;  % Relative deviation of the cost f.
   if cost<cost0
      X0=X;cost0=cost;              % Updating X0 (best multi-sine spectra)
      iter0=iter;
      relax=relax/2;                % Lowering Levenberg-Marquardt factor
   else
      X=X0;cost=cost0;              % Restoring the best results
      x=f2t(X,N); 
      if ~isempty(H), Y=(H.*X)*W(2)/W(1); y=f2t(Y,N); end
      relax=relax*10;               % Augmenting Levenberg-Marquardt factor
   end
   
   fprintf('P=%g/CF=%g/ITER=%g',p,crestfactor(X,Nmax,Fe),iter)
   fprintf('/NORM=%g/ITER0=%g/NORM0=%g',cost,iter0,cost0)
   fprintf('/LM=%g\n',relax)
end

X=X0/W(1);  % Unscaling the result
