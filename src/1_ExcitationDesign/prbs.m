function [bitseries,X,freq,nextstnum]=prbs(fs,log2N,bitno,startnum)
%PRBS - pseudo-random binary sequence, PRBS.
%   [bitseries,X,freq,nextstnum]=prbs(log2N,bitno,startnum,fs)
% bitseries     : generated bit series (values +1,-1).
%                 Column vector, length: bitno, default length: 2^log2N-1
% nextstnum     : startnum to be used if continued sequence generation
% log2N         : bit length of the shift register (int 2<=log2N<=30)
% bitno         : number of bits to be generated (length of bitseries)
% startnum      : digital equivalent of shift register start value (opt)
% Algorithm     : Maximum Length Binary sequence (MLBS)
% Author        : Thomas Beauduin, KULeuven, 2014
%%%%%
if nargin<3, bitno=2^log2N-1; end
if nargin<4, startnum=2^log2N-1; end

startnum=round(rem(abs(startnum-1),2^log2N-1))+1;
stn=startnum; reg=zeros(log2N,1);
for i=1:log2N
  reg(i)=rem(stn,2); 
  stn=(stn-reg(i))/2; 
end

zind=find(reg==0); reg(zind)=-1*ones(length(zind),1);
if     log2N==2, multind=[1,2];
elseif log2N==3, multind=[2,3];
elseif log2N==4, multind=[3,4];
elseif log2N==5, multind=[3,5];
elseif log2N==6, multind=[5,6];
elseif log2N==7, multind=[4,7];
elseif log2N==8, multind=[4,5,6,8];
elseif log2N==9, multind=[5,9];
elseif log2N==10, multind=[7,10];
elseif log2N==11, multind=[9,11];
elseif log2N==12, multind=[6,8,11,12];
elseif log2N==13, multind=[9,10,12,13];
elseif log2N==14, multind=[4,8,13,14];
elseif log2N==15, multind=[14,15];
elseif log2N==16, multind=[4,13,15,16];
elseif log2N==17, multind=[14,17];
elseif log2N==18, multind=[11,18];
elseif log2N==19, multind=[14,17,18,19];
elseif log2N==20, multind=[17,20];
elseif log2N==21, multind=[19,21];
elseif log2N==22, multind=[21,22];
elseif log2N==23, multind=[18,23];
elseif log2N==24, multind=[17,22,23,24];
elseif log2N==25, multind=[22,25];
elseif log2N==26, multind=[20,24,25,26];
elseif log2N==27, multind=[22,25,26,27];
elseif log2N==28, multind=[25,28];
elseif log2N==29, multind=[27,29];
elseif log2N==30, multind=[7,28,29,30];
end

bitseries=zeros(bitno,1);
for i=1:bitno
  bitseries(i)=-prod(reg(multind));
  reg=[bitseries(i);reg(1:log2N-1)];
end
nextstnum=sum((reg/2+1/2).*(2.0.^(0:log2N-1)'));

X = t2f(bitseries,length(bitseries));
freq = fs*(0:1:length(bitseries)/2-1)'/length(bitseries);
end