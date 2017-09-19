function [Y] = dbm(X)
%DBM - Decibel magnitude

Y = 20*log10(abs(X));


end

