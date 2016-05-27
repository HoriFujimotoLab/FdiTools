function X = msinl2p(X,nrofs,itp,Fa,W,H)
%MSINL2P - Lp-norm optimization of multisine phases.
%
% X         : Matrix with input spectra (1=DC)
% nrofs     : Number of samples per period in time domain
% Fa        : additional harmonic numbers for snow multisine
% W         : Weighting row vector, e.g., W=[1=same_cr, 2=same_a];
% H         : Transfer function matrix for input-output opti
% iterno    : Max. number of iterations
% relvar    : Max. relative variation of L2p-norm
% itp       : initial guess of phases (optional)
% X         : Matrices with optimal input spectra
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
if nargin<6, H=[]; end              % in-only minimization
if nargin<5, W=0; end               % in-out opti weighting
if nargin<4, Fa=[]; end             % Additional snowing freq
W=W(:)'; Fa=Fa(:);                  % vectorization
p = 2;                              % initial l2p-norm value

Fe = find(X~=0);                    % effective spectral lines
Ft = [Fe; Fa];                      % total spectral lines
nrofe = length(Fe);                 % number of effective harm
nrofa = length(Fa);                 % number of additional harm
nroft = length(Ft);                 % number of considered harm
kmax = max(Ft)-1;                   % Largest harmonic number
if length(X(:,1)) <= kmax, X(kmax+1,1) = 0; end

% Parameter Initialization
if strcmp(itp,'r'), 
    tex = 'Random';                 % limit opti to keep random 
    relvar = 1e-6;                  % max relative deviation
    iterno = 10;                    % max number of iterations
else
    tex = 'Schroed';                % maximize opti for good cf
    relvar = 1e-10;                 % max relative deviation
    iterno = 20;                    % max number of iterations
end

% Gradually increase the p-norm
fprintf('\n Iterative calculation: min. L2P-norm \n');
while p < 500

% Identity between L2p-summation and integration criterions
Ndummy = 2*p*kmax+1;      
N = 2; while Ndummy >= N, N = N*2; end
N = min([N, nrofs]);

% Input (and output) crest factor compression
if ~isempty(H)
   H = H(1:kmax+1); Y = H.*X;
   if W==0, W = [effval(X,Fe),effval(Y,Fe)].^(-1); end % Same crest factors
   if W==1, W = [effval(X,Fe),effval(X,Fe)].^(-1); end % Same extreme value
   X = X*W(1); Y = Y*W(2);
   x = f2t(X,N);               % Input signals
   y = f2t(Y,N);               % Output signals
else
   W = 1/effval(X,Fe);
   X = X*W;                    % Scaling: max|x(t)|~1
   x = f2t(X,N);               % Input signals
end

if isempty(H),  cost = lpnorm(x(:),2*p);
else            cost = lpnorm([x(:);y(:)],2*p);
end

% Optimization initial values
relax = 0.01;                           % Levenberg-Marquardt
iter = 0;                               % Iteration counter
relerror = inf;                         % Relative cost deviation
X0=X; cost0=cost; iter0=0;              % Init best results
indexmat = (Ft-1)*ones(size(Ft))';        
indexmins = indexmat - indexmat' + N;   % Init of (k-l)-matrix
indexplus = indexmat + indexmat' + N;   % init of (k+l)-matrix
phkphlmins = indexmins(:);              % Matrices to vector
phkphlplus = indexplus(:);
if ~isempty(Fa)
   phkalmins = indexmins(:,nrofe+1:nroft);
   phkalplus = indexplus(:,nrofe+1:nroft);
   akalmins = indexmins(nrofe+1:nroft,nrofe+1:nroft);
   akalplus = indexplus(nrofe+1:nroft,nrofe+1:nroft);
   phkalmins = phkalmins(:); akalmins = akalmins(:);
   phkalplus = phkalplus(:); akalplus = akalplus(:);
end

% MAIN L2P-norm minimization loop
while (iter<iterno)&&(relerror>relvar)
   iter = iter+1;
   X2p2 = fft(x.^(2*p-2));
   X2p2 = [flipud(conj(X2p2(2:N)));X2p2];
   
   % Calculation of JphJph and JphE
   Qmins = conj(X(Ft)*X(Ft)');
   Qplus = conj(X(Ft)*X(Ft).');
   Qmins(:) = Qmins(:).*X2p2(phkphlmins);
   Qplus(:) = Qplus(:).*X2p2(phkphlplus);
   JphJph = p*real(Qmins - Qplus);
   JphE = sum(imag(Qmins + Qplus)')';
   
   if ~isempty(Fa)
       % Calculation of JphJa and JaE
       Qmins = conj(X(Ft)*exp(1i*angle(X(Fa)))');
       Qplus = conj(X(Ft)*exp(1i*angle(X(Fa))).');
       Qmins(:) = Qmins(:).*X2p2(phkalmins);
       Qplus(:) = Qplus(:).*X2p2(phkalplus);
       JphJa = p*imag(Qmins+Qplus);
       JaE = sum(real(Qmins+Qplus))';
      
       % Calculation of JaJa
       Qmins = conj(exp(1i*angle(X(Fa)))*exp(1i*angle(X(Fa)))');
       Qplus = conj(exp(1i*angle(X(Fa)))*exp(1i*angle(X(Fa))).');
       Qmins(:) = Qmins(:).*X2p2(akalmins);
       Qplus(:) = Qplus(:).*X2p2(akalplus);
       JaJa = p*real(Qmins + Qplus);
   end

   if ~isempty(H)
       Y = (H.*X)*W(2)/W(1);
       y = f2t(Y,N);
       Y2p2 = fft(y.^(2*p-2));
       Y2p2 = [flipud(conj(Y2p2(2:N)));Y2p2];

       % Calculation of JphJph and JphE
       Qmins=conj(Y(Ft)*Y(Ft)');
       Qplus=conj(Y(Ft)*Y(Ft).');
       Qmins(:)=Qmins(:).*Y2p2(phkphlmins);
       Qplus(:)=Qplus(:).*Y2p2(phkphlplus);
       JphJph=JphJph+p*real(Qmins-Qplus);
       JphE=JphE+sum(imag(Qmins+Qplus)')';
      
       if ~isempty(Fa)
           % Calculation of JphJa and JaE
           Qmins = conj(Y(Ft)*exp(1i*angle(Y(Fa)))');
           Qplus=conj(Y(Ft)*exp(1i*angle(Y(Fa))).');
           Qmins(:)=Qmins(:).*Y2p2(phkalmins);
           Qplus(:)=Qplus(:).*Y2p2(phkalplus);
           JphJa=JphJa+p*imag(Qmins+Qplus);
           JaE=JaE+sum(real(Qmins+Qplus))';

           % Calculation of JaJa
           Qmins = conj(exp(1i*angle(Y(Fa)))*exp(1i*angle(Y(Fa)))');
           Qplus = conj(exp(1i*angle(Y(Fa)))*exp(1i*angle(Y(Fa))).');
           Qmins(:) = Qmins(:).*Y2p2(akalmins);
           Qplus(:) = Qplus(:).*Y2p2(akalplus);
           JaJa=JaJa + p*real(Qmins + Qplus);
       end
   end

   % Least Squares matrices
   if isempty(Fa),  A = JphJph; b = JphE;
   else             A = [JphJph,JphJa; JphJa',JaJa]; b = [JphE;JaE];
   end     
   A = A + relax * diag(diag(A)+max(diag(A))*eps);

   % Levenberg-Marquardt Algorithm (updating)
   Delta = -A\b;
   X(Fe) = abs(X(Fe)).*exp(1i*(angle(X(Fe))+Delta(1:nrofe)));
   X(Fa) = (abs(X(Fa))+Delta(nroft+1:nroft+nrofa))...
           .*exp(1i*(angle(X(Fa))+Delta(nrofe+1:nroft)));
   x = f2t(X,N);
   if ~isempty(H),
       Y = (H.*X)*W(2)/W(1); 
       y = f2t(Y,N);
   end
   
   if isempty(H),   cost = lpnorm(x(:),2*p);            % L2p-norm calc
   else             cost = lpnorm([x(:);y(:)],2*p);
   end
   
   relerror = abs(cost-cost0)/cost0;        % Relative cost deviation
   if cost < cost0
      relax = relax/2;                      % Lowering LM factor
      X0 = X; cost0 = cost; iter0 = iter;   % Updating best results
   else
      relax = relax*10;                     % Augmenting LM factor 
      X = X0; cost = cost0;                 % Restoring best results
      x = f2t(X,N); 
      if ~isempty(H), 
          Y = (H.*X)*W(2)/W(1);
          y = f2t(Y,N); 
      end
   end
   CF = lpnorm(x,inf)./effval(X,Fe);
   fprintf('P=%g/iter=%g: CF=%g, err=%g\n',p,iter,CF,relerror)
end
X = X0/W(1);            % Unscaling the result
p = ceil(p*2);          % Increasing the p-norm

end
fprintf('\n CREST FACTOR = %g \n\n',CF)

end
