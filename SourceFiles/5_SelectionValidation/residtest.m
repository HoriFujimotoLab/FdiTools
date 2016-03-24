function [Lags,Corr,CB,Fraction,Tag] = residtest(freq,FRFs,FRF,sCR,nrofp)
% RESID - Identification Residuals for validation
%
% freq      : freq vector of measurement
% FRFs      : FRF measurement
% FRF       : FRF of fitted model
% sCR       : Cramer-Rao underbound of measurement
% norfp     : nr of periods in measurment
% Lags      : lags for auto_corr ploting
% Corr      : Auto correlation
% CB        : Confindence Bounds
% Fraction  : fraction under confidence bounde
% Author    : Thomas Beauduin, KULeuven, 2014
%
% see AUTO-CORRELATION figure book, fig.11-8

nroff = length(freq);
Res=zeros(1,nroff);
Fraction = zeros(1,2);
Lags = (-nroff+1:nroff-1);
Tag = fieldnames(FRF);
nrofm = length(Tag);

% Confidence Bounds
scale0 = (nrofp-2)/(nrofp-1);
scale = (nrofp-5/3)/(nrofp-11/12);
Conf_scale0=scale0*((nrofp-1)^(3/2)/(nrofp-2)/(nrofp-3)^(1/2));
Conf_scale=scale*(nrofp-1)/(nrofp-2);

% std auto-correlation
p50 = sqrt(-log(1-0.5));
p95 = sqrt(-log(1-0.95));
ConfScale = Conf_scale*ones(size(Lags));
ConfScale(nroff) = Conf_scale0;
Conf_Bound = repmat(ConfScale./(nroff-abs(Lags)).^(0.5),[2,1]);
Conf_Bound(1,:) = p50*Conf_Bound(1,:);
Conf_Bound(2,:) = p95*Conf_Bound(2,:);
CB=Conf_Bound;

% Auto-correlation
FRF_c = struct2cell(FRF);
Select = (1:2*nroff-1);
Select(nroff) = [];
for m = 1:nrofm
    FRFm = FRF_c{m};
    for f = 1:nroff
        Res(f) = (FRFs(f)-FRFm(f))./sCR(f).^0.5;    % Residuals
    end
    Auto_Corr = xcorr(Res,'unbiased');  % ./nroff auto-correlation
    TheScale = scale*ones(size(Lags));
    TheScale(nroff) = scale0;
    Auto_Corr = Auto_Corr .* TheScale;
    Corr(m,:)=squeeze(abs(Auto_Corr));         % abs-corr (plot)
    
    Fraction(m,1) = length(find(abs(Auto_Corr(Select)) - ...
                          Conf_Bound(1,Select) > 0)) / (2*nroff-2)*100;
    Fraction(m,2) = length(find(abs(Auto_Corr(Select)) - ...
                          Conf_Bound(2,Select) > 0)) / (2*nroff-2)*100;
end

[Fraction(:,1),index]=sort(Fraction(:,1),'descend');
Fraction(:,2)=Fraction(index,2);
Corr=Corr(index,:);
Tag = Tag(index);
end

