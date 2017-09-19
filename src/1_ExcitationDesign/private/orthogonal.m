function T = orthogonal(nrofi)
%ORTHOGONAL - Transformation matrix for orthogonal multisines.
%   T = orthogonal(nrofi)
% nrofi : Number of inputs (odd inputs only).
% T     : Create nrofi orthogonal phase multisines
% Author: Thomas Beauduin, KULeuven PMA, 2014
%%%%%

for p=1:nrofi
    for q=1:nrofi
        T(p,q)=nrofi^(-1/2)*exp(1i*2*pi*(p-1)*(q-1)/nrofi);
    end
end
%T = T./T(1,1);

end