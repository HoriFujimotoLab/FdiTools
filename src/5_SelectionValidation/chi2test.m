function [err,var,tag] = chi2test(X,Y,freq,FRF,sCR,SYS)
% CHI2TEST - Identification Chi-Squares test (MIMO).
%   [err,var,tag] = chi2test(X,Y,freq,FRF,sCR,SYS)
% X,Y,freq  : Input & output frequency domain data
% FRF       : Measurement Frequncy Response Matrix
% sCR       : Cramer-Rao underbound of measurement
% SYS       : structured estimated models for comparison
% err,tag   : absolute error between tag model and data
% var       : absolute 95% chi^2 confidence bounds
% Author    : Thomas Beauduin, KULeuven, PMA division, 2014
%%%%%
nrofi = size(X,2);                  % number of inputs
nrofo = size(Y,2);                  % number of outputs
nroff = length(freq(:));            % number of frequency lines
nrofh = nrofi*nrofo;                % number of transfer functions
nrofm = length(struct2cell(SYS));   % number of system models

% Calculation of 95% chi2 confidence bounds
model_c = struct2cell(SYS);
N_alfa = 10.5966;
var = zeros(nroff,nrofh);
for h=1:nrofh
    var(:,h) = N_alfa*sCR(:,h)/2;
end

% Calculation of residual modelling error
tag = cell(nrofm,nrofh);
err = zeros(nroff,nrofm,nrofh);
srt = zeros(nroff,nrofm,nrofh);
cnt = zeros(nrofm,nrofh);
for h=1:nrofh
    i = ceil(h/nrofo); o = h-(i-1)*nrofo;
    for m=1:nrofm
        [Bn,An] = tfdata(model_c{m}(o,i),'v');
        FRFsys = squeeze(freqresp(tf(Bn,An),freq*2*pi));
        err(:,m,h) = abs(FRF(:,h)-FRFsys).^2;
        for f=1:nroff
            if err(f,m,h) > var(f,h)
                cnt(m,h) = cnt(m,h)+1;
                srt(f,m,h) = err(f,m,h) - var(f,h);
            end
        end
    end
    [~,index] = sort(mean(srt(:,:,h),1).*cnt(:,h)','descend');
    err(:,:,h) = err(:,index,h);
    tag(:,h) = fieldnames(SYS);
    tag(:,h) = tag(index,h);
end
