function [Hbtls,Hgtls] = btlsfdi(varargin)
%BTLSFDI - Bootstrapped Total Least Squares Estimation (MIMO).
% structured input
%   [Hbtls,Hgtls] = btlsfdi(Pest,n,mh,ml,iter,relax,max_err)
% Pest        : Estimated model structure (frd) obtained by time2frf_ml.m
% FdiTools classical input
%   [Hbtls,Hgtls] = btlsfdi(X,Y,freq,n,mh,ml,sY2,sX2,cXY,iter,relax,max_err,fs)
% X,Y,freq    : Input & output frequency domain data
% sX2,sY2     : variance of X & Y frequency domain data
% cXY         : Covariance between X & Y frequency domain data
% n,mh,ml     : Order of the denominator/nominator polynomials
% relax       : Relaxation factor: 0 =< relax =< 1 (r=0: GTLS, r=1: full BTLS)
% max_iter    : Maximum number of iterations (stop criterion)
% max_err     : Maximum model relative error (stop criterion)
% cORd, fs    : Continuous or discrete time model identification
% Hbtls,Hgtls : BTLS/GTLS iterative & initial estimation solution
% Author      : Thomas Beauduin, KULeuven, PMA division, 2014
%               Wataru Ohnishi, The University of Tokyo, 2019
%%%%%

if length(varargin) < 9 % structured input
    Pest = varargin{1};
    n = varargin{2};
    M_mh = varargin{3};
    M_ml = varargin{4};
    relax = varargin{5};
    max_iter = varargin{6};
    max_err = varargin{7};
    cORd = varargin{8};

    X = Pest.UserData.X;
    Y = Pest.UserData.Y;
    freq = Pest.freq;
    sX2 = Pest.UserData.sX2;
    sY2 = Pest.UserData.sY2;
    cXY = Pest.UserData.cXY;
    if iscell(Pest.UserData.ms), fs = Pest.UserData.ms{1}.harm.fs;
    else fs = Pest.UserData.ms(1).harm.fs; end
else % FdiTools classical input
    X = varargin{1};
    Y = varargin{2};
    freq = varargin{3};
    n = varargin{4};
    M_mh = varargin{5};
    M_ml = varargin{6};
    sY2 = varargin{7};
    sX2 = varargin{8};
    cXY = varargin{9};
    relax = varargin{10};
    max_iter = varargin{11};
    max_err = varargin{12};
    cORd = varargin{13};
    fs = varargin{14};
end


nrofi = size(X,2);                      % number of inputs
nrofo = size(Y,2);                      % number of outputs
nrofh = nrofi*nrofo;                    % number of transfer functions
nroff = length(freq(:));                % number of frequency lines
nrofb = sum(M_mh-M_ml)+nrofh;           % number of numerator coefficients
nrofp = nrofb+n;                        % number of estimated parameters
M_mh=M_mh'; M_ml=M_ml';                 % vectorize numerator sizes
M_mh = M_mh(:); M_ml = M_ml(:);

% Calculation of initial values for iterative process
fprintf(' \n Initial calculation: GTLS solution \n')
[Hgtls,waxis] = gtlsfdi(X,Y,freq,n,M_mh,M_ml,sX2,sY2,cXY,cORd,fs);

% Calculation of iterative parameter estimation
fprintf('\n Iterative calculation: BTLS solution \n');
[Bb,Ab] = hm2ba(Hgtls);                 % starting values choice
iter = 0;                               % interation number
err = max_err + 1;                      % model relative error
Xg0 = [1;ba2theta(Bb,Ab,n,M_mh,M_ml)];

while (iter<=max_iter)&&(err>max_err)
    iter = iter+1;
    Den = polyval(Ab,waxis);
    Num = zeros(nroff,nrofh);
    for h=1:nrofh, Num(:,h) = polyval(Bb(h,:),waxis); end
    
    % Calculation of weighted Jacobian (WreJre(Z))
    Xh = zeros(nroff,nrofh); 
    Yh = zeros(nroff,nrofh);
    for h=1:nrofh
        i = ceil(h/nrofo); o = h-(i-1)*nrofo;
        SEr = sqrt((sX2(:,i).*(abs(Num(:,h)).^2)...
                  + sY2(:,o).*(abs(Den).^2)...
                  - 2*real(cXY(:,h).*Den.*conj(Num(:,h)))).^relax);
        Xh(:,h) = X(:,i)./SEr; Yh(:,h) = Y(:,o)./SEr;
    end
    EX = kron(ones(nrofh*nroff,1),(n:-1:0));
    W = kron(ones(nrofh,n+1),waxis);
    P = (W.^EX).*kron(ones(1,n+1),Yh(:));
    
    index = 1;
    Q = zeros(nrofh*nroff,nrofb);
    for h=1:nrofh
        EX = kron(ones(nroff,1),(M_mh(h):-1:M_ml(h)));
        W = kron(ones(1,M_mh(h)-M_ml(h)+1),waxis);
        U = (W.^EX).*kron(ones(1,M_mh(h)-M_ml(h)+1),Xh(:,h));
        Q(nroff*(h-1)+1:nroff*h,index:index+M_mh(h)-M_ml(h)) = U;
        index = index + M_mh(h)-M_ml(h)+1;
    end
    J = [real(P) -real(Q); imag(P) -imag(Q)];

    % Calculation of column covariance matrix (Cwj^1/2)
    C = zeros(nrofh*(nrofp+1),nrofp+1);
    for h=1:nrofh
        i=ceil(h/nrofo); o=h-(i-1)*nrofo;
        MytMy = zeros(n+1,n+1);
        MxtMx = zeros(M_mh(h)-M_ml(h)+1,M_mh(h)-M_ml(h)+1);
        MytMx = zeros(n+1,M_mh(h)-M_ml(h)+1);
        SEr2 = ( sX2(:,i).*(abs(Num(:,h)).^2)...
              + sY2(:,o).*(abs(Den).^2)...
              - 2*real(cXY(:,h).*Den.*conj(Num(:,h)))).^relax;
        for p=n:-1:0
            for q=n:-1:0
                MytMy(n-p+1,n-q+1) = ...
                2*real(((-1)^q)*sum(waxis.^(p+q).*(sY2(:,o)./SEr2)));
            end
        end
        for p=M_mh(h):-1:M_ml(h)
            for q=M_mh(h):-1:M_ml(h)
                MxtMx(M_mh(h)-p+1,M_mh(h)-q+1) = ...
                2*real(((-1)^q)*sum(waxis.^(p+q).*(sX2(:,i)./SEr2)));
            end
        end
        for p=n:-1:0
            for q=M_mh(h):-1:M_ml(h)
                MytMx(n-p+1,M_mh(h)-q+1) = ...
                -2*real(((-1)^q)*sum(waxis.^(p+q).*(cXY(:,h)./SEr2)));
            end
        end
        MtM=[MytMy MytMx ; MytMx' MxtMx];
        cols = (1:n+1+M_mh(h)-M_ml(h)+1);
        rows = (h-1)*(nrofp+1)+1:(h-1)*(nrofp+1)+(n+1+M_mh(h)-M_ml(h)+1);
        try chol(A)
            C(rows,cols) = chol(MtM);
        catch 
            C(rows,cols) = [chol(MytMy) MytMx ; MytMx' MxtMx];
            disp('Matrix is not positive definite')
        end   
    end

    % Calculation of generalized right singular vector (Xg)
    Xg = qsvd(J,C);
    Xg = inv(Xg');
    Xg = Xg(:,nrofp+1)/Xg(1,nrofp+1);
    [Bb,Ab] = theta2ba(Xg(2:end),n,M_mh,M_ml);
    err = max(abs((Xg0 - Xg)./Xg));
    Xg0 = Xg;
    cost = btlsfdi_res(Bb,Ab,freq,X,Y,sX2,sY2,cXY,relax,waxis);
    
    fprintf('Iter %g: index = %g, cost = %g, rel.err = %g\n',...
    iter,iter,cost,err)  
end
Hbtls = ba2hm(Bb,Ab,nrofi,nrofo);

end