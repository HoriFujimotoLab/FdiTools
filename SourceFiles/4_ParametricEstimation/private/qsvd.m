function [X,Ua,Ub,Sa,Sb]=qsvd(A,B)
% QSVD - Generalized Singular Value Decomposition.
%   [X,Ua,Ub,S1,S2]=qsvd(A,B)
% A, B  : input matrices A[n x m] - B[p x m]
% Ua,Ub : Unitary Matrices Ua[n x m] - Ub[p x m] 
% Sa,Sb : Singular Diagonal Matrix: Sa[diag(s1,..,sm)]
% X     : Generalized right singular vectors
%         A = Ua*Sa*X' , B = Ub*Sb*X'
% Author: Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
% Important Notes:
% 1. This function can only be applied to matrices
%    A and B that have more rows than columns.
[m,n1]=size(A); [p,n2]=size(B);
    if      n1 == n2, n=n1;
    else    error('equal number of columns needed')
    end
rab = n; ra = min(n,m); rb = min(n,p);

[U,S,V] = svd([A;B]);
U1 = U(1:m,1:rab);
U2 = U(m+1:m+p,1:rab);

[Ua,S1,V1] = svd(U1);
S2 = sqrt(eye(rab)-S1'*S1);
Ub = (U2*V1*pinv(S2));

Sa=[ eye(rab-rb),zeros(rab-rb,rb);
     zeros(ra+rb-rab,rab-rb),S1(1:ra+rb-rab,1:ra+rb-rab),zeros(ra+rb-rab,n-ra);
	 zeros(m-ra,n)];
Sb=[ zeros(p-rb,n);
     zeros(ra+rb-rab,rab-rb),S2(1:ra+rb-rab,1:ra+rb-rab),zeros(ra+rb-rab,n-ra)
	 zeros(rab-ra,ra),eye(rab-ra) zeros(rab-ra,n-rab)];
X=(V*S(1:n,1:n)*V1);

end