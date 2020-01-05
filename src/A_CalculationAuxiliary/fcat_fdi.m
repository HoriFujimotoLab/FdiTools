function sysout = fcat_fdi_test(varargin)
% FdiTools version of fcat
% Wataru Ohnishi, The University of Tokyo, 2019
%%%%

N = length(varargin);
idx_frd = false(N,1);
for k = 1:N
    tmp = varargin{k};
    tmp2 = whos('tmp');
    if strcmp(tmp2.class,'frd'), idx_frd(k) = true; end
end
if sum(~idx_frd) > 1
    error('error in option');
elseif sum(~idx_frd) == 1
    option = varargin(~idx_frd);
    option = option{1};
else
    option.noise = 'sGhat';
end
frddata = varargin(idx_frd);

Nfrd = length(frddata);
name = fieldnames(frddata{1}.UserData);
idx_ms = strcmp(name,'ms'); % field 'ms' is not frequency dependent
idx_x = strcmp(name,'x'); % field 'x' is not frequency dependent
idx_y = strcmp(name,'y'); % field 'y' is not frequency dependent
idx_nofreq = or(idx_ms,idx_x);
idx_nofreq = or(idx_nofreq,idx_y);

sysout = frddata{1};
for k = 2:Nfrd
    temp = sysout;
    sysout = fcat_fdi_resp(sysout,frddata{k},option);
    %     sysout = fcat_fdi_UserData(name(~idx_nofreq),sysout,temp,frddata{k});
end

if isfield(frddata{1}.UserData,'ms') % multisine
    Nfrd2 = Nfrd;
    for k = 1:Nfrd
        Nfrd2 = Nfrd2 + length(frddata{k}.UserData.ms) - 1;
    end
    sysout.UserData.ms = cell(1,Nfrd2);
    kk = 1;
    for k = 1:Nfrd
        if iscell(frddata{k}.UserData.ms)
            sysout.UserData.ms(kk:kk+length(frddata{k}.UserData.ms)-1) = ...
                frddata{k}.UserData.ms(kk:kk+length(frddata{k}.UserData.ms)-1);
            kk = kk + length(frddata{k}.UserData.ms);
        else
            sysout.UserData.ms{kk} = frddata{k}.UserData.ms;
            kk = kk + 1;
        end
    end
end

if exist('Nfrd2')
    sysout.UserData.x = cell(1,Nfrd2);
    sysout.UserData.y = cell(1,Nfrd2);
    kk = 1;
    for k = 1:Nfrd
        if isfield(frddata{k}.UserData,'x') % time domain data exists
            if iscell(frddata{k}.UserData.x)
                sysout.UserData.x(kk:kk+length(frddata{k}.UserData.x)-1) = ...
                    frddata{k}.UserData.x(kk:kk+length(frddata{k}.UserData.x)-1);
                sysout.UserData.y(kk:kk+length(frddata{k}.UserData.y)-1) = ...
                    frddata{k}.UserData.y(kk:kk+length(frddata{k}.UserData.y)-1);
                kk = kk + length(frddata{k}.UserData.x) -1;
            else
                sysout.UserData.x{kk} = frddata{k}.UserData.x;
                sysout.UserData.y{kk} = frddata{k}.UserData.y;
            end
        end
        kk = kk + 1;
    end
end


end

function sysout = fcat_fdi_resp(sys1,sys2,option)
if length(sys1.freq) < length(sys2.freq)
    temp = sys1;
    sys1 = sys2;
    sys2 = temp; clear temp
end

freq2 = sys2.freq;
for k = 1:length(freq2)
    k1 = find(sys1.freq==freq2(k));
    if k1
        k2 = sys2.freq==freq2(k);
        noise1 = getfield(sys1.UserData,option.noise);
        noise2 = getfield(sys2.UserData,option.noise);
        if noise1(k1) < noise2(k2)
            sys2 = fdel_fdi(sys2,freq2(k),freq2(k));
        else
            sys1 = fdel_fdi(sys1,freq2(k),freq2(k));
        end
    end
end

sysout = fcat_fdi_UserData(sys1,sys2);

end

function sysout = fcat_fdi_UserData(varargin)

N = length(varargin);

sysout = varargin{1};
name = fieldnames(sysout.UserData);
idx_ms = strcmp(name,'ms'); % field 'ms' is not frequency dependent
idx_x = strcmp(name,'x'); % field 'x' is not frequency dependent
idx_y = strcmp(name,'y'); % field 'y' is not frequency dependent
idx_nofreq = or(idx_ms,idx_x);
idx_nofreq = or(idx_nofreq,idx_y);
vars = name(~idx_nofreq);

for k = 2:N
    temp = sysout;
    sysout = fcat(sysout,varargin{k});
    sysout = fcat_fdi_UserData2(vars,sysout,temp,varargin{k});
end

% if isfield(varargin{1}.UserData,'ms') % multisine
%     sysout.UserData.ms = cell(1,N);
%     for k = 1:N
%         sysout.UserData.ms{k} = varargin{k}.UserData.ms;
%     end
% end
% if isfield(varargin{1}.UserData,'x') % time domain data
%     sysout.UserData.x = cell(1,N);
%     sysout.UserData.y = cell(1,N);
%     for k = 1:N
%         sysout.UserData.x{k} = varargin{k}.UserData.x;
%         sysout.UserData.y{k} = varargin{k}.UserData.y;
%     end
% end

end

function sysout = fcat_fdi_UserData2(sname,sysin,sys1,sys2)
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


% function sysout = fcat_fdi_UserData(sname,sysin,sys1,sys2)
% N = length(sname);
% for k = 1:N
%     temp1 = getfield(sys1.UserData,sname{k});
%     temp2 = getfield(sys2.UserData,sname{k});
%     freqmerge = [sys1.freq;sys2.freq;];
%     fieldmerge = [temp1;temp2;];
%     [~,I] = sort(freqmerge);
%     fieldmerge = fieldmerge(I,:);
%     %     sysin.UserData = rmfield(sysin.UserData,sname{k});
%     sysin.UserData = setfield(sysin.UserData,sname{k},fieldmerge);
% end
% sysout = sysin;
% end

