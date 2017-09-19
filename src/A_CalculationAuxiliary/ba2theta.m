function y = ba2theta(Bn,An,n,M_mh,M_ml)
%BA2THETA - rational polynomial to parameter vector (MIMO).
%    y = ba2theta(Bn,An,n,M_mh,M_ml)
% y         : estimation parameter vector (theta_est)
% n,M_mh/ml : transfer function polynomial orders
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%
nrofh = length(M_ml(:));        

y(1:n,1) = An(2:n+1);

index = n+1;
for i=1:nrofh
   yy=Bn(i,n+1-M_mh(i):n+1-M_ml(i));
   y(index:index + M_mh(i)-M_ml(i))=yy;
   index = index + M_mh(i)-M_ml(i) + 1;
end

end

