function Pest = fdicohere(Pest)
% calculate coherence for fdi object created by time2frf_ml
% Wataru Ohnishi, the University of Tokyo, 2019
%%%%

nrofs = Pest.UserData.ms.nrofs;
fs = Pest.UserData.ms.harm.fs;
[cxy,freq] = mscohere(Pest.UserData.x,Pest.UserData.y,...
    rectwin(nrofs),0,nrofs,fs);

if isfield(Pest.UserData.ms,'ex')
    Pest.UserData.cxy = cxy(Pest.UserData.ms.ex,:);
else
    [~,kmin] = min(abs(freq-Pest.UserData.ms.harm.fl));
    [~,kmax] = min(abs(freq-Pest.UserData.ms.harm.fh));
    Pest.UserData.cxy = cxy(kmin:kmax,:);
end

end
