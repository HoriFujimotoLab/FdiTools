function [cost,interval,tag] = costtest(X,Y,freq,sX2,sY2,cXY,SYS,relax,nrofp)
% COSTTEST - Modelling residual cost validation test.
%   [cost,interval,tag] = costtest(X,Y,freq,sX2,sY2,cXY,SYS,relax,nrofp)
% X,Y,freq  : Input & output frequency domain data
% sX2,sY2   : variance of X & Y frequency domain data
% cXY       : Covariance between X & Y frequency domain data
% SYS       : Structure containing the estimated models
% norfp     : number of periods per measurment
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
nrofi = size(X,2);                  % number of inputs
nrofo = size(Y,2);                  % number of outputs
nrofh = nrofi*nrofo;                % number of transfer functionss
nrofm = length(struct2cell(SYS));   % number of system models

% Calculation of interval
model_c = struct2cell(SYS);
n = length(tfdata(model_c{1}(:,1),'v'))-1;
Vnoise =((nrofp-1)/(nrofp-2))*(length(freq)-n);
interval = [Vnoise-2*Vnoise^(1/2) Vnoise+2*Vnoise^(1/2)];

% Calculation of residual cost
tag = cell(nrofm,nrofh);
cost = zeros(nrofm,nrofh);
for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    for m=1:nrofm
        [Bn,An] = tfdata(model_c{m}(o,i),'v');
        cost(m,h) = fdicost(Bn,An,freq,X(:,i),Y(:,o),...
                            sX2(:,i),sY2(:,o),cXY(:,h),relax);
    end
    tag(:,h) = fieldnames(SYS);
end
[~,index] = sort(sum(cost,2),'descend');
cost = cost(index,:);
tag = tag(index,:);

end

% Note:
% 1. this is for estimator selection
%    create additional function for model set/order selection
%    input: data + matrix of parameters + requirested estimator
