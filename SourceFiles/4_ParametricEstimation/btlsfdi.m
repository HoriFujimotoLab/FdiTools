function [Bb,Ab,Bg,Ag] = btlsfdi(X,Y,freq,n,mh,ml,sY2,sX2,cXY,iter,relax,max_dev)
%BTLSFDI Bootstrapped Total Least Squares Estimation (SISO).
%
% X,Y       : input & output values of the FRF
% freq      : frequency vector
% sX2, sY2  : variance of X & Y frequency domain noise
% cXY       : covariance between X & Y frequency domain noise
% n         : order of the denominator polynomial
% mh, ml    : high & low order of the numerator polynomial
% relax     : relaxation factor : 0 =< relax =< 1 (r=0: GTLS, r=1: full BTLS)
% iter      : number of BLTS iterations
% Bb,Ab     : BTLS solution
% Bg,Ag     : GTLS solution
% max_dev   : maximum model relative deviation (stop criterion)
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%%
% GTLS initial guess
[Bg,Ag,xqsvd] = gtlsfdi(Y,X,freq,n,mh,ml,sY2,sX2,cXY);

freq=freq(:);
w = 1i*freq*2*pi;
nA = (n:-1:0);
nB = (mh:-1:ml);
Bb = Bg;
Ab = Ag;
xqsvd_o = xqsvd;

%start of iterative BTLS
tel=1; mod_change = max_dev+1;
while (tel<=iter)&&(mod_change>max_dev)
    Num = polyval(Bb,w);
    Den = polyval(Ab,w);
    scal2=(sX2.*(abs(Num).^2)+sY2.*(abs(Den).^2)...
           -2*real(cXY.*Den.*conj(Num))).^relax;
    scal = sqrt(scal2);
    Ys=Y./scal; Xs=X./scal;
    sX2s = sX2./scal2;
    sY2s = sY2./scal2;
    cXYs = cXY./scal2;

    A = (kron(w,ones(size(nA)))...
        .^kron(ones(size(w)),nA)).*kron(Ys,ones(size(nA)));
    A = [A -(kron(w,ones(size(nB)))...
        .^kron(ones(size(w)),nB)).*kron(Xs,ones(size(nB)))];
    A = [real(A);imag(A)];

    MytMy = zeros(n+1,n+1);
    MxtMx = zeros(mh-ml+1,mh-ml+1);
    MytMx = zeros(n+1,mh-ml+1);

    for (i=n:-1:0)
        for (j=n:-1:0)
            MytMy(n-i+1,n-j+1) = real(((-1)^j)*2*sum(w.^(i+j).*sY2s));
        end;
    end;

    for (i=mh:-1:ml)
        for (j=mh:-1:ml)
            MxtMx(mh-i+1,mh-j+1) = real(((-1)^j)*2*sum(w.^(i+j).*sX2s));
        end;
    end;

    for (i=n:-1:0)
        for (j=mh:-1:ml)
            MytMx(n-i+1,mh-j+1) =  -real(((-1)^j)*2*sum(w.^(i+j).*cXYs));
        end;
    end;

    MtM=[MytMy MytMx ; MytMx' MxtMx];
    M = chol(MtM);
  
    %qsvd solution
    [Ua,Un,Sa,Sn,xqsvd]=qsvd(A,M);
    nt=n+1+mh-ml+1;
    xqsvd=inv(xqsvd');
    xqsvd=xqsvd(:,nt)/xqsvd(1,nt);
    Ab = xqsvd(1:n+1)';
    Bb = [zeros(1,n-mh) xqsvd(n+2:nt)' zeros(1,ml)];
    mod_change=max(abs((xqsvd_o-xqsvd)./xqsvd));
    xqsvd_o=xqsvd;

    Num = polyval(Bb,w);
    Den = polyval(Ab,w);
    T = abs(Num.*X - Den.*Y).^2;
    scal2_n = ( sX2.*(abs(Num).^2) + sY2.*(abs(Den).^2) ...
              - 2*real(cXY.*Den.*conj(Num)) );
    cost = sum(T./scal2)/sum(scal2_n./scal2);
    fprintf('BTLS-iteration %g of %g; cost = %g; max. model change %g \n'...
            ,tel,iter,cost,mod_change)
    tel = tel+1;
end
end