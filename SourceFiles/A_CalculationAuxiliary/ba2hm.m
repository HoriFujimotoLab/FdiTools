function [Hm] = ba2hm(Bn,An,nrofi,nrofo)
%BA2HM - parameter matrices to transfer function array.
%
% Bn,An     : nominator and denominator coefficient matrices
% nrofi/o   : number of input & outputs of the system
% SYS       : Transfer function matrix [nrofi x nrofo]
% Author    : Thomas Beauduin, KULeuven, PMA, 2014
%%%%%
nrofh = nrofi*nrofo;             % number of transfer functions

for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    [z,p,k] = tf2zpk(Bn(h,:),An);
    Hm(o,i) = zpk(z,p,k);
end

end

