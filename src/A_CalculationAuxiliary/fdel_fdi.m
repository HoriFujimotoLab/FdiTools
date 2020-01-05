function sysout = fdel_fdi(sys, fmin, fmax, vars)
% FdiTools version of fdel
% Wataru Ohnishi, The University of Tokyo, 2019
%%%%

if nargin < 4
    name = fieldnames(sys.UserData);
    idx_ms = strcmp(name,'ms'); % field 'ms' is not frequency dependent
    idx_x = strcmp(name,'x'); % field 'x' is not frequency dependent
    idx_y = strcmp(name,'y'); % field 'y' is not frequency dependent
    idx_nofreq = or(idx_ms,idx_x);
    idx_nofreq = or(idx_nofreq,idx_y);
    vars = name(~idx_nofreq);
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

