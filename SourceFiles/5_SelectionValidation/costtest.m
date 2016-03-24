function [cost,interval,tag] = costtest(SYS,X,Y,freq,sX2,sY2,cXY,cORd,fs,relax,nrofp)
% COSTTEST - Modelling residual cost validation test.
%
% freq      : freq vector of measurement
% FRFm      : FRF measurement
% FRF       : FRF of fitted model
% sCR       : Cramer-Rao underbound of measurement
% norfp     : nr of periods in measurment
% Lags      : lags for auto_corr ploting
% Corr      : Auto correlation
% CB        : Confindence Bounds
% Fraction  : fraction under confidence bounde
% Author    : Thomas Beauduin, KULeuven, 2014
%
model_n = fieldnames(SYS);
model_c = struct2cell(SYS);
B=tfdata(model_c{1},'v');
n=length(B)-1;

Vnoise =((nrofp-1)/(nrofp-2))*(length(freq)-n);
interval=[Vnoise-2*Vnoise^(1/2) Vnoise+2*Vnoise^(1/2)];
for k=1:length(model_c)
    [B,A]=tfdata(model_c{k},'v');
    cost(k)=btlsfdi_res(B,A,freq,X,Y,sX2,sY2,cXY,cORd,fs,relax);
end
[cost,index]=sort(cost,'descend');
tag = model_n(index);

end


% Note:
% this is for estimator selection
% create additional function for model set/order selection
% input: data + matrix of parameters + requirested estimator
