function y = mkn(x,cval )
%MKN replaces 0's in a mask with NaN
%   RJO wrote to simpify making nanmasks
%    written 4/29/14

y = double(x);

if nargin == 1
    y(x==0)=NaN;
    
elseif cval > 0
        y(x==0) = NaN + NaN*1i;
else
        y(x==0)=NaN;
end

end

