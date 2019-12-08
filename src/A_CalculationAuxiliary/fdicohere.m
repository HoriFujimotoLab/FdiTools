function Pest = fdicohere(Pest)
% calculate coherence for fdi object created by time2frf_ml
% Wataru Ohnishi, the University of Tokyo, 2019
%%%%

nrofs = length(Pest.UserData.ms.x);
fs = Pest.UserData.ms.harm.fs;
[cxy,~] = mscohere(Pest.UserData.x,Pest.UserData.y,...
    rectwin(nrofs),0,nrofs,fs);
Pest.UserData.cxy = cxy(Pest.UserData.ms.ex,:);

end
