function [Y] = phs(varargin)
%PHS - phase of frequency domain data.
%   [Y] = phs(X)
% X,Y   : input/output imaginary vectors
% glitch: trigger level for glitch removal
% method: shift removal and glitch removal:
%         shift : up-down jump patterns in wrapped data
%         glitch: unnaturaly strong diff in data
% author: Thomas Beauduin, University of Tokyo, 2016

X = varargin{1};
Y = angle(X)*180/pi;

% 1. SHIFT REMOVAL
% remove wrapped data up-down shifts for visibility
shift = 160; j=0; 
first = 0; last=0;
for i=2:length(Y)
    jump_i=0; jump_j=0; end_j=0;
    if abs(Y(i)-Y(i-1))>shift
        if j~=0
            if Y(i)>Y(i-1), jump_i=+1;  % jump up
            else            jump_i=-1;
            end
        else                bgn_j=+1;   % begin of array
        end

        for j=i:length(Y)-1
            if abs(Y(j)-Y(j+1))>shift, break; end
        end
        
        if j~=length(Y)-1
            if Y(j)>Y(j+1), jump_j=-1;  % jump down
            else            jump_j=+1;
            end
        else                end_j=+1;   % end of array
        end
        
        if bgn_j == 1 && first == 0
            if Y(i)<Y(i-1), Y(1:i-1)=Y(1:i-1)-360;%abs(Y(i)-Y(i-1));
            else            jump_i=+1;
            end
            first = 1;
        end
        if jump_i == -jump_j
            Yij = mean(Y(i:j));
            Y(i:j)=Y(i:j)-360;%abs(Yij-Y(i-1));%mean([abs(Yij-Y(i-1)),abs(Yij-Y(j+1))]);
        end
        if end_j == 1 && last == 0
           if nargin > 1, Y(i:j+1)=Y(i:j+1)-360;%abs(Y(i)-Y(i-1)); 
           end
           last = 1;
        end
    end
end

% 2. GLITCH REMOVAL
% remove cleaned data up-down glitch for visibility
if nargin > 2
    glitch = varargin{end};
    for i=3:length(Y)-2
        if ( abs(Y(i)-Y(i-1))>glitch || abs(Y(i)-Y(i-2))>glitch ) && ...
           ( abs(Y(i)-Y(i+1))>glitch || abs(Y(i)-Y(i+2))>glitch )
            Y(i)=mean([Y(i-3),Y(i-2),Y(i-1),Y(i+1),Y(i+2),Y(i+3)]);
        end
    end
end


% % 2. clean phase glitches
% glitch = 90;
% for i=3:length(Y)-2
%     if ( abs(Y(i)-Y(i-1))>glitch || abs(Y(i)-Y(i-2))>glitch ) && ...
%        ( abs(Y(i)-Y(i+1))>glitch || abs(Y(i)-Y(i+2))>glitch )
%         %Y(i)=mean([Y(i-3),Y(i-3),Y(i-2),Y(i-2),Y(i-1),Y(i-1),...
%         %           Y(i+1),Y(i+2),Y(i+3)]);
%         Y(i)=Y(i-3);
%     end
% end

% if abs(Y(i)-Y(i-1)) > abs(Y(j)-Y(j+1))
%     Y(i:j)=Y(i:j)-abs(Y(i)-Y(i-1));
% else
%     Y(i:j)=Y(i:j)-abs(Y(j)-Y(j+1));
% end


end

