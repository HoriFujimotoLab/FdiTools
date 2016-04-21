function [sys_c] = sysclean(sys_i)
%SYSCLEAN - Clean identified model based on requirements.
%
% sys_i     : identified system model
% sys_c     : cleaned system model
% Author    : Thomas Beauduin, University of Tokyo, 2016
%%%%%
[zp,pp,kp] = zpkdata(sys_i,'v');
nrofh = length(pp);                 % Number of tf's
nrofp = length(pp{1});              % Number of poles

index = 0;
for h=1:nrofh       % number of tf
    for p=1:nrofp   % number of pp
        if real(pp{h}(p)) > 0, pp{h}(p) = 0; end
        if abs(pp{h}(p)) == 0, index = 1;    end
    end
    if index == 0, 
        pp{h}(end) = 0; 
        pp{h}(end-1) = real(pp{h}(end-1)); 
    end
    

    nrofz = length(zp{h});
    for z=1:nrofz   % number of zp
        if real(zp{h}(z)) > 0, zp{h}(z) = -real(zp{h}(z))+1i*imag(zp{h}(z)); end
    end
end
sys_c = zpk(zp,pp,kp);

end

