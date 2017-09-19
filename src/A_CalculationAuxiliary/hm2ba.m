function [Bn,An] = hm2ba(Hm)
%HM2BA - parameter matrices to transfer function array.
%
% Bn,An     : nominator and denominator coefficient matrices
% nrofi/o   : number of input & outputs of the system
% SYS       : Transfer function matrix [nrofi x nrofo]
% Author    : Thomas Beauduin, KULeuven, PMA, 2014
%%%%%
nrofi = size(Hm,2);             % number of inputs
nrofo = size(Hm,1);             % number of outputs 
nrofh = nrofi*nrofo;            % number of tf's

for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    [Bn(h,:),An] = tfdata(Hm(o,i),'v');
end

end

