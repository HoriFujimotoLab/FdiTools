function [fqlog,fidx] = lin2qlog(flin,fr)
%LIN2QLOG - Linear to quasi-logarithmic frequency grid.
%
% fqlog     : Quasi-logarithmic frequency lines
% fidx      : Frequency lines index in linear grid
% Author    : Thomas Beauduin, KULeuven, 2014
%%%%%
last = 1; fidx(1) = 1;          % frequency index's
fqlog(1) = flin(last);          % quasi-log freq grid

while last < length(flin)
    [~,next] = min(abs(flin(last)*fr - flin(last+1:length(flin))));
    last = last + next;
    
    % if step is acceptably large
    if flin(last) > fqlog(length(fqlog))*sqrt(fr),
        fidx = [fidx; last];
        fqlog = [fqlog; flin(last)];
    end
end

end

