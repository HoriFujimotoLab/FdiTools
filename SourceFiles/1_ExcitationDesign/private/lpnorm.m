function Lp = lpnorm(A,p)
% LPNORM - Calculate the vector Lp-norm.
%   Lp = lpnorm(A,p)
% A         : matrix of column vectors
% p         : number of the norm
% Author    : Thomas Beauduin, KULeuven, 2014
%%%%%
if nargin==1, p=2; end

[nRow,nCol] = size(A);
for i=1:nCol, Lp(i) = norm(A(:,i),p); end
if p~=inf, Lp = Lp/(nRow^(1/p)); end
