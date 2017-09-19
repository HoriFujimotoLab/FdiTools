function [A,B,C,D,fs]=ssfdi(FRF,FRF_W,freq,n_in,n_out,n_row,n_col,order)
% (NOTE: WORK IN PROGRESS)
% Frequency domain identification of a discrete time MIMO state space model
% method of T. McKelvey (subspace approach)
%
% FRF   : matrix with FRF [H11 H21 H31 ... H12 H22 H32 ...]  
% FRF_W : weighting matrix for FRF
% freq  : frequency axis maximum frequency value becomes Nyquist freq of model
% n_in  : number of inputs
% n_out : number of outputs
% n_row : number of block rows in Hankel matrix
% n_col : number of block columns in Hankel matrix
% order : order of resulting model (OPTIONAL parameter)
% A,B,C,D : discrete state space model 
% fs      : sampling frequency of model, is equal to max(freq)*2
%%%%%

freq=freq(:);
nrofmea = length(freq);
fs=2*max(freq);
FRF=[FRF;conj(FRF(nrofmea-1:-1:2,:))];
FRF(nrofmea,:)=real(FRF(nrofmea,:));
nr = nrofmea;
nrofmea = 2*nrofmea-2;
omega = 2*pi*(0:1:nrofmea-1)'/nrofmea;
ifrf=real(ifft(FRF));

H=[];
for (i=1:n_in);
  q=ifrf(:,n_out*(i-1)+1:n_out*i)';q=q(:);
  H=[H q];
end;

  
Hqr=[];
for (i=1:n_col)
  Hqr=[Hqr H(n_out*i+1:n_out*(i+n_row),:)];
end;

[UU,SS,VV]=svd(Hqr);

if (nargin==8)
  n = order;
else
  semilogy(diag(SS),'*');grid
  n=input('the system order equals : ');
end;
    

U1=UU(:,1:n);
V1=VV(:,1:n);
S1=SS(1:n,1:n);
Li = U1*sqrt(S1);

C= Li(1:n_out,:);
A = pinv(Li(1:n_out*(n_row-1),:))*Li(n_out+1:n_out*n_row,:);

H=[];
for (i=1:n_in);
  q=conj(FRF(:,n_out*(i-1)+1:n_out*i))';q=q(:);
  H=[H q];
end;
Wght=[];
for (i=1:n_in);
  q=conj(FRF_W(:,n_out*(i-1)+1:n_out*i))';q=q(:);
  Wght=[Wght q];
end;


S=[];
for (i=1:nr)
  CA = C*inv(exp(j*omega(i))*eye(n,n)-A); 
  S = [S; [eye(n_out,n_out) CA]];
end;
S=[real(S);imag(S)];


D=[];
B=[];
for (i=1:n_in)
  W=kron([Wght(1:n_out*nr,i);Wght(1:n_out*nr,i)],ones(1,n_out+n));
  DB = pinv(W.*S)*(W(:,1).*[real(H(1:n_out*nr,i)); imag(H(1:n_out*nr,i))]);
  D = [D DB(1:n_out,:)];
  B = [B DB(n_out+1:n_out+n,:)];  
end;    
