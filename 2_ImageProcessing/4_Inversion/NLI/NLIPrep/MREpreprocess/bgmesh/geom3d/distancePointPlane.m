function d = distancePointPlane(point, plane)
%DISTANCEPOINTPLANE : compute euclidean distance betwen 3D point and plane
%
%   D = distancePointPlane(POINT, PLANE) return the distance between point
%   POINT and the plane PLANE, given as :
%   POINT : [x0 y0 z0]
%   PLANE : [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
%   D     : scalar  
%   
%   ---------
%
% Songbai Ji (6/29/2006). Allow one plane, many pointz; many planes one
% point; or N planes and N points configuration in the input.
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY


if size(plane, 1) == 1; %one plane possible many points
    plane = repmat(plane, size(point,1), 1);
elseif size(point,1) == 1; % one point and many planes
    point = repmat(point, size(plane,1),1);
elseif (size(plane,1) ~= size(point,1)) ; % N planes and M pointa, not allowed for now.
    error('input size not correct, either one/many plane and many/one line, or same # of planes and lines!');
end


% normalized plane normal
n = normalize(cross(plane(:,4:6), plane(:, 7:9), 2));


% Uses hessian form, ie : N.p = d
% In this case, d can be found as : -N.p0, when N is normalized
d = -dot(n, plane(:,1:3)-point(:,1:3), 2);

