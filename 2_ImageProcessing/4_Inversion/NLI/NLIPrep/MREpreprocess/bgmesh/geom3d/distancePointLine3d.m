function [d] = distancePointLine3d(point, line)
%distancePointLine3d : compute euclidean distance betwen 3D point and line
%
%   D = distancePointLine3d(POINT, LINE) return the distance between point
%   POINT and the plane LINE, given as :
%   POINT : [x0 y0 z0]
%   PLANE : [x0 y0 z0 dx dy dz]
%   D     : scalar  
%   
%   ---------
% 
% Modifiedy by Songbai Ji (6/26/2006).  Now also allow one point, many
% lines; many points one line; or N points and N lines configuration in the
% input.
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 23/05/2005.
%

%   HISTORY

if size(point, 1) == 1; %one plane possible many lines
    point = repmat(point, size(line,1), 1);
elseif size(line,1) == 1; % one line and many planes
    line = repmat(line, size(point,1),1);
elseif (size(point,1) ~= size(line,1)) ; % N planes and M lines, not allowed for now.
    error('input size not correct, either one/many point and many/one line, or same # of points and lines!');
end
    
% cf. Mathworld (distance point line 3d)  for formula
% d = norm(cross(line(:,4:6), (line(:,1:3)-point)))./norm(line(:,4:6));
d = normByRow(cross(line(:,4:6), (line(:,1:3)-point)))./normByRow(line(:,4:6));

function nrmByRow = normByRow(m)
% helper function to compute norm for each row of the matrix m.
%  Right now uses for loop, maybe there is a vecterized way???
% 
% Songbai Ji (6/26/2006).
rows = size(m,1);
nrmByRow = zeros(rows,1);
for i = 1 : rows
    nrmByRow(i) = norm(m(i,:));
end
