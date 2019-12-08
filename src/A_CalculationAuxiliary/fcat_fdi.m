function sysout = fcat_fdi(varargin)
% FdiTools version of fcat
% Wataru Ohnishi, The University of Tokyo, 2019
%%%%

N = length(varargin);

sysout = varargin{1};
for k = 2:N
    temp = sysout;
    sysout = fcat(sysout,varargin{k});
    sysout = fcat_fdi_UserData({'X','Y','FRFn','sX2','sY2','cXY','sCR'},sysout,temp,varargin{k});
end

sysout.UserData.ms = varargin{1}.UserData.ms;
for k = 2:N
    sysout.UserData.ms = [sysout.UserData.ms;varargin{k}.UserData.ms;];
end

end

function sysout = fcat_fdi_UserData(sname,sysin,sys1,sys2)
N = length(sname);
for k = 1:N
    temp1 = getfield(sys1.UserData,sname{k});
    temp2 = getfield(sys2.UserData,sname{k});
    freqmerge = [sys1.freq;sys2.freq;];
    fieldmerge = [temp1;temp2;];
    [~,I] = sort(freqmerge);
    fieldmerge = fieldmerge(I,:);
%     sysin.UserData = rmfield(sysin.UserData,sname{k});
    sysin.UserData = setfield(sysin.UserData,sname{k},fieldmerge);
end
sysout = sysin;
end

