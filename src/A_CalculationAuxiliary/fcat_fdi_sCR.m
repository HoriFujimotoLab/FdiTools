function sysout = fcat_fdi_sCR(sys1,sys2)
% FdiTools version of fcat considering sCR
% Wataru Ohnishi, The University of Tokyo, 2019
%%%%

if length(sys1.freq) < length(sys2.freq)
    temp = sys1; 
    sys1 = sys2;
    sys2 = temp; clear temp
end

freq = sys2.freq;
for k = 1:length(freq)
    k1 = find(sys1.freq==freq(k));
   if k1
       k2 = find(sys2.freq==freq(k));
        if sys1.UserData.sCR(k1) < sys2.UserData.sCR(k2)
            sys2 = fdel_fdi(sys2,freq(k),freq(k));
        else
            sys1 = fdel_fdi(sys1,freq(k),freq(k));
        end
   end
end

sysout = fcat_fdi(sys1,sys2);

end

