function Lp = lpnorm(A,p)
% LPNORM - Calculate the vector Lp-norm.
%
% A         : matrix of column vectors
% p         : number of the norm
% Lp        : Lp norm
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
% see also MSINL2PI, MSINCLIP

if nargin==1, p=2; end

[nRow,nCol] = size(A);
for i=1:nCol
    Lp(i) = norm(A(:,i),p);
end
if p~=inf
    Lp = Lp/(nRow^(1/p));
end
