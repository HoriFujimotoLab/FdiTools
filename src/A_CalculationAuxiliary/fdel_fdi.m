function sysout = fdel_fdi(sys, fmin, fmax, vars)
% FdiTools version of fdel
% Wataru Ohnishi, The University of Tokyo, 2019
%%%%

if nargin < 4
    if isfield(sys.UserData,'sCR'), flag_ms = true; else, flag_ms = false; end
    if isfield(sys.UserData,'cxy'), flag_ch = true; else flag_ch = false; end
    if flag_ms && ~flag_ch
        vars = {'X','Y','FRFn','sX2','sY2','cXY','sCR'};
    elseif flag_ms && flag_ch
        vars = {'X','Y','FRFn','sX2','sY2','cXY','sCR','cxy'};
    else
        vars = {'cxy'};
    end
end

[~,kmin] = min(abs(sys.freq - fmin));
[~,kmax] = min(abs(sys.freq - fmax));

sysout = fdel(sys,sys.freq(kmin:kmax));
sysout = fdel_fdi_UserData(sysout,vars, kmin,kmax);

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

