function sysout = fcat_fdi(varargin)
% FdiTools version of fcat
% Wataru Ohnishi, The University of Tokyo, 2019
%%%%

N = length(varargin);

sysout = varargin{1};
if isfield(sysout.UserData,'sCR'), flag_ms = true; else, flag_ms = false; end
if isfield(sysout.UserData,'cxy'), flag_ch = true; else flag_ch = false; end
if flag_ms && ~flag_ch
    vars = {'X','Y','FRFn','sX2','sY2','cXY','sCR'};
elseif flag_ms && flag_ch
    vars = {'X','Y','FRFn','sX2','sY2','cXY','sCR','cxy'};
else
    vars = {'cxy'};
end

for k = 2:N
    temp = sysout;
    sysout = fcat(sysout,varargin{k});
    sysout = fcat_fdi_UserData(vars,sysout,temp,varargin{k});
end

if isfield(varargin{1}.UserData,'ms') % multisine
    sysout.UserData.ms = cell(1,N);
    for k = 1:N
        sysout.UserData.ms{k} = varargin{k}.UserData.ms;
    end
end
if isfield(varargin{1}.UserData,'x') % time domain data
    sysout.UserData.x = cell(1,N);
    sysout.UserData.y = cell(1,N);
    for k = 1:N
        sysout.UserData.x{k} = varargin{k}.UserData.x;
        sysout.UserData.y{k} = varargin{k}.UserData.y;
    end
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

