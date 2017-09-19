function [Bn,An] = theta2ba(y,n,M_mh,M_ml)
%THETA2BA - parameter vector to rational polynomial (MIMO).
%   [Bn,An] = theta2ba(y,n,M_mh,M_ml)
% y         : estimation parameter vector (theta_est)
% n,M_mh/ml : transfer function polynomial orders
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%
nrofh = length(M_ml(:));

An = [1 y(1:n)'];
Bn = zeros(nrofh,n+1);
index = n+1;
for i=1:nrofh
   yy = y(index:index + M_mh(i)-M_ml(i))';
   Bn(i,n+1-M_mh(i):n+1-M_ml(i))=yy;
   index = index + M_mh(i)-M_ml(i) + 1;
end

end