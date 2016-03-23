function [Ua,Ub,S1,S2,X]=qsvd(A,B)
% QSVD - Quadratic Singular Value Decomposition: M=U*S*X'
%
% A, B  : input matrices
% Ua,Ub : Unitary Matrices 
% S1,S2 : Recangular Diagonal Matrix
% X     : Transposed Unitary Matrix
% Author: Thomas Beauduin, KULeuven, 2014
%
% Important Notes:
% 1. This function can only be applied to matrices
%    A and B that have more rows than columns.
% 2. Ub is not square (as in the exact definition).  
%    The zeros in Sb and the columns in Ub corresponding to these
%    zeros are missing (they are hardly ever needed any way).

[m,n1]=size(A);
[p,n2]=size(B);

if n1 ~= n2, error('equal number of columns needed')
else         n=n1;
end

[U S V]=svd([A;B]);

rab=n;
ra=min(n,m);
rb=min(n,p);
U1=U(1:m,1:rab);
U2=U(m+1:m+p,1:rab);
[Ua S1 V1]=svd(U1);
S2=sqrt(eye(rab)-S1'*S1);
%Ub=[zeros(p,p-rab) U2*V1*pinv(S2)];
Ub=[U2*V1*pinv(S2)];

Sa=[	eye(rab-rb) zeros(rab-rb,rb);
zeros(ra+rb-rab,rab-rb) S1(1:ra+rb-rab,1:ra+rb-rab) zeros(ra+rb-rab,n-ra);
	zeros(m-ra,n)];
Sb=[	zeros(p-rb,n);
zeros(ra+rb-rab,rab-rb) S2(1:ra+rb-rab,1:ra+rb-rab) zeros(ra+rb-rab,n-ra)
	zeros(rab-ra,ra) eye(rab-ra) zeros(rab-ra,n-rab)];
X=[V*S(1:n,1:n)*V1];

end
