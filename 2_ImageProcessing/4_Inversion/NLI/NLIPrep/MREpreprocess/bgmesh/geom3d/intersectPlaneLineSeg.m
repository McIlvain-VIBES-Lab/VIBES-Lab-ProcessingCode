function points = intersectPlaneLineSeg(plane,x,y,z)
% function points = intersectPlaneLineSeg(plane,x,y,z);
% 
% Program that finds the intersections of a plane and line segments denoted
% by 3D points [x y z].  
% 
% NOTE: need to check if the hashing-unique scheme below actually works.
% 
% Songbai Ji (01/15/2007).

%% check input
x=x(:);y=y(:);z=z(:);
if length(unique([length(x) length(y) length(z)]))~=1
	error('input x, y and z must be of same size!');
end

%% create the lines and find intersections with plane
pnt1 = [x(1:end-1) y(1:end-1) z(1:end-1)];
pnt2 = [x(2:end) y(2:end) z(2:end)];

lines = createLine3d(pnt1, pnt2);

points = intersectPlaneLine(plane, lines);

%% check if point is within the line segment, if so, found
cen = [mean([pnt1(:,1) pnt2(:,1)],2), mean([pnt1(:,2) pnt2(:,2)],2), ...
	mean([pnt1(:,3) pnt2(:,3)],2)];

tmp = pnt2-pnt1;
segLength = sqrt(sum(tmp.*tmp,2));

tmp = points - cen;
p2cen = sqrt(sum(tmp.*tmp,2));

% if intersection to center distance less or equal to segment length, then
% good:
points = points( p2cen<=segLength, : );

%% but possibly the line segment end points are counted > once.  If this
% happend, we only choose one:
% create hash function
hash = sin(points(:,1)) + sin(points(:,2)) + sin(points(:,3)) ;

[hash,I]=unique(hash);
points = points (I,:);

