function [confid,var,tag] = chi2test(FRFs,FRF,sCR)
% CHI2TEST - Chi-Squares test for identification validation.
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

model_n = fieldnames(FRF);
model_c = struct2cell(FRF);

for k=1:length(model_c)
    confid(k,:)=(abs(FRFs-model_c{k}).^2)./1;
end
N_alfa=10.5966;                % 95% chi2 confidence bounds
var=N_alfa*sCR/2;

tag = model_n;

cost = mean(confid,2);
[~,index]=sort(cost,'descend');
tag = model_n(index);
confid = confid(index,:);