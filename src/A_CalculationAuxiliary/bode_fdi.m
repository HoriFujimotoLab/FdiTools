function hfig = bode_fdi(data,noise,option)

if nargin < 3
    option.pmin = -180;
    option.pmax = 180;
end

N = length(data); % number of data
freq = logspace(0,3,400);
for k = 1:N
    try freq = data{k}.freq; catch, data{k} = frd(data{k},freq,'FrequencyUnit','Hz'); end
end

hfig = figure;
subplot(2,1,1);
for k = 1:N
    h = semilogx(data{k}.frequency,mag2db(abs(squeeze(data{k}.ResponseData)))); hold on;
end

subplot(2,1,2);
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
end

subplot(2,1,1);
h = semilogx(noise(:,1),mag2db(abs(noise(:,2)))); hold on;

end