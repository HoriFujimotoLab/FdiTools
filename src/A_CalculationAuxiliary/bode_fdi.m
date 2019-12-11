function hfig = bode_fdi(data,noise,option)
% number of plot
if ~iscell(data)
    data = {data};
end
N = length(data); % number of data

Nplot = 3;
cohFlag = true;
for k = 1:N
    if isfield(data{k}.UserData,'sCR'), Nplot = 2; cohFlag = false; end
end

if nargin < 2
    if cohFlag, noise = 'cxy'; else, noise = 'FRFn'; end
end

if nargin < 3
    option = [];
end
if ~isfield(option,'pmin'), option.pmin = -180; option.pmax = 180; end

freq = logspace(0,3,400);
for k = 1:N
    try freq = data{k}.freq; catch, data{k} = frd(data{k},freq,'FrequencyUnit','Hz'); end
end
hfig = figure;
subplot(Nplot,1,1);
for k = 1:N
    h = semilogx(data{k}.frequency,mag2db(abs(squeeze(data{k}.ResponseData)))); hold on;
end
if isfield(option,'title'), title(option.title); end
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

if cohFlag
    subplot(Nplot,1,3);
    for k = 1:N
        h = semilogx(data{k}.frequency,data{k}.UserData.cxy); hold on;
    end
    ylabel('Coherence [-]');
    xlabel(['Frequency [',data{1}.FrequencyUnit,']']);
else
    xlabel(['Frequency [',data{1}.FrequencyUnit,']']);
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
    if ~cohFlag
        if N == 1
            legend('FRFs',noise);
        else
            strLeg = cell(N*2,1);
            for k = 1:N
                strLeg{k} = sprintf('FRFs%d',k); strLeg{k+N} = sprintf('FRFn%d',k); 
            end
            legend(strLeg);
        end
    end
end


end
