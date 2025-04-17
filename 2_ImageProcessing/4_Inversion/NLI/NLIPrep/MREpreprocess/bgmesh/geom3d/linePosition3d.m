function d = linePosition3d(point, line)
%LINEPOSITION3D return position of a 3D point on a 3D line
%
%   L = LINEPOSITION3D(POINT, LINE)
%   compute position of point POINT on the line LINE, relative to origin
%   point and direction vector of the line.
%   LINE has the form [x0 y0 dx dy],
%   POINT has the form [x y], and is assumed to belong to line.
%
%   L = LINEPOSITION(POINT, LINES)
%   if LINES is an array of NL lines, return NL positions, corresponding to
%   each line.
%
%   L = LINEPOSITION(POINTS, LINE)
%   if POINTS is an array of NP points, return NP positions, corresponding
%   to each point.
%
%   L = LINEPOSITION(POINTS, LINES)
%   if POINTS is an array of NP points and LINES is an array of NL lines,
%   return an array of [NP NL] position, corresponding to each couple
%   point-line.
%
%   see createLine for more details on line representation.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%

%   HISTORY :

if size(point, 1)== 1 && size(line, 1)>1
    point = repmat(point, [size(line, 1) 1]);
end

if size(point, 1)>1 && size(line, 1)==1
    line = repmat(line, [size(point, 1) 1]);
end


dp = point - line(:,1:3);
dl = line(:,4:6);

d = dot(dp, dl, 2)./dot(dl, dl, 2);
