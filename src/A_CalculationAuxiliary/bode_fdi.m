function hfig = bode_fdi(data,noise,option)
% number of plot
if isfield(data.UserData,'cxy'), Nplot = 3; else, Nplot = 2; end

if nargin < 2
    if isfield(data.UserData,'cxy'), noise = 'cxy'; else, noise = 'sCR'; end
end

if nargin < 3
    option.pmin = -180;
    option.pmax = 180;
end

if ~iscell(data)
    data = {data};
end

N = length(data); % number of data
freq = logspace(0,3,400);
for k = 1:N
    try freq = data{k}.freq; catch, data{k} = frd(data{k},freq,'FrequencyUnit','Hz'); end
end

hfig = figure;
subplot(Nplot,1,1);
for k = 1:N
    h = semilogx(data{k}.frequency,mag2db(abs(squeeze(data{k}.ResponseData)))); hold on;
end
ylabel('Magnitude [dB]');

subplot(Nplot,1,2);
for k = 1:N
    phasedeg = rad2deg(angle(squeeze(data{k}.ResponseData)));
    for kk = 1:length(phasedeg)
        while phasedeg(kk) > option.pmax
            phasedeg(kk) = phasedeg(kk) - 360;
        end
        while phasedeg(kk) < option.pmin
            phasedeg(kk) = phasedeg(kk) + 360;
        end
    end
    h = semilogx(data{k}.frequency,phasedeg); hold on;
    yticks(option.pmin:90:option.pmax);
    ylim([option.pmin,option.pmax]);
end
ylabel('Phase [deg]');

if isfield(data{k}.UserData,'cxy')
    subplot(Nplot,1,3);
    for k = 1:N
        h = semilogx(data{k}.frequency,data{k}.UserData.cxy); hold on;
    end
    ylabel('Coherence [-]');
else
    subplot(2,1,1);
    if ischar(noise)
        for k = 1:N
            if isfield(data{k}.UserData,noise)
                h = semilogx(data{k}.freq,mag2db(abs(getfield(data{k}.UserData,noise)))); hold on;
            end
        end
    else
        h = semilogx(noise(:,1),mag2db(abs(noise(:,2)))); hold on;
    end
end

xlabel(['Frequency [',data{1}.FrequencyUnit,']']);

end
