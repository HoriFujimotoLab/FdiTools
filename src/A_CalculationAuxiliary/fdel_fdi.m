function sysout = fdel_fdi(sys, fmin, fmax)
% FdiTools version of fdel
% Wataru Ohnishi, The University of Tokyo, 2019
%%%%

[~,kmin] = min(abs(sys.freq - fmin));
[~,kmax] = min(abs(sys.freq - fmax));

sysout = fdel(sys,sys.freq(kmin:kmax));
sysout = fdel_fdi_UserData(sysout,{'X','Y','FRFn','sX2','sY2','cXY','sCR'}, kmin,kmax);

end

function out = fdel_fdi_UserData(in,sname, kmin,kmax)
N = length(sname);
for k = 1:N
    temp = getfield(in.UserData,sname{k});
    in.UserData = rmfield(in.UserData,sname{k});
    in.UserData = setfield(in.UserData,sname{k},[temp(1:kmin-1,:);temp(kmax+1:end,:);]);
end
out = in;
end

