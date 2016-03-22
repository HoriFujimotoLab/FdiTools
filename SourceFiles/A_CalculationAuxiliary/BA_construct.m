function [Bn,An]=BA_construct(y,n,M_mh,M_ml)
%BA_CONSTRUCT build transfer function form.
% 
% y         : periodic broad-band excitation measurement
% n         : transfer function numerator order
% M_mh/_ml  : lowest & highest frequencies in band
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
%%%%
nino=length(M_ml(:));
maxB = max(n+1,max(n-min(M_ml))+1);

% storing the solution in matrices An and Bn
An = [1 y(1:n)'];
Bn = zeros(nino,maxB);
index_count = n+1;
for (i=1:nino)
   yy = y(index_count:index_count + M_mh(i)-M_ml(i))';
   Bn(i,n+1-M_mh(i):n+1-M_ml(i))=yy;
 
   index_count = index_count + M_mh(i)-M_ml(i) + 1;
end;

% VERSION 2:
%   maxB = max(n+1,max(max(M_mh))+1);
%   Bn(i,maxB-M_mh(i):maxB-M_ml(i))=yy;