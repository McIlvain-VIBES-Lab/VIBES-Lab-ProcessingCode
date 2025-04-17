function [x, y, I] = circlesort(sx, sy, start)
% [x, y, I] = circlesort(sx, sy [, start]). 
% 
% Sort the data points to be in a counter-clock wise order. Input as column
% vectors, start is the starting angle relative to horizonal line (in
% radians).  If the point cloud is 3D, use the projected data on one of the
% main axis planes.
%
%
% NOTE: only good for data convex to the centroid of the point cloud.
%
% Songbai Ji (6/25/2006): also return a 3rd output: I, such that x=sx(I);
%  
% See also closestsort.m;

if nargin <3;
    start = 0; 
else
    
end;
mx = mean(sx); my = mean(sy);
[th, r] = cart2pol (sx-mx, sy-my);  %-pi<th<pi

th = th - start; % srat from the named line
Ineg = find(th<0);
th(Ineg) = th(Ineg) + 2*pi;  % now: 0 - 2*pi

[th_s, I] = sort(th);

x= sx(I); y = sy(I);














% older version:
% x=[]; y=[];
% if length(sx) ~= length(sy)
%     disp('length of inputs must be the same');
%     return;
% end
% [r, c] = size(sx); if r == 1 sx = sx'; end;
% [r, c] = size(sy); if r == 1 sy = sy'; end;
% 
% cenx= sum(sx)/length(sx);
% ceny = sum(sy)/length(sy);
% 
% I = find (sy-ceny>=0);
% upx = sx(I); upy = sy(I);
% I = find(sy-ceny<0);
% downx = sx(I); downy = sy(I);
% upx2 = upx - cenx; upy2 = upy - ceny;
% downx2 = downx -cenx; downy2 = downy - ceny;
% %first the upper plane:
% diag = sqrt(upx2.^2 + upy2.^2);
% cosine = upx2./diag;
% [Y, I]=sort(cosine);
% for i = 1: length(I)
%     I2(i) = I(length(I) - i +1);
% end
% 
% %now the lower plane:
% diag = sqrt(downx2.^2 + downy2.^2);
% cosine = downx2./diag;
% [Y, I]=sort(cosine);
% x = [upx(I2); downx(I)];
% y = [upy(I2); downy(I)];
