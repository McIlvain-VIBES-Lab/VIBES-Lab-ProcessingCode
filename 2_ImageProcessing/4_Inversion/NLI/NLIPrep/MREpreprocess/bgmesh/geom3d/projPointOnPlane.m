function point = projPointOnPlane(point, plane)
%PROJPOINTONPLANE : return the projection of a point on a plane
%
%   PT2 = PROJECTEDPOINT(PT1, PLANE).
%   Compute the (orthogonal) projection of point PT1 onto the PLANE.
%   
%   Function works also for multiple points and planes. In this case, it
%   returns multiple points.
%   Point PT1 is a [N*3] array, and PLANE is a [N*9] array (see createPlane
%   for details). Result PT2 is a [N*3] array, containing coordinates of
%   orthogonal projections of PT1 onto planes PLANE.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%
%   Corrected by Songbai Ji, 1/11/2008.

n = planeNormal(plane);
if size(point,1)==1 && size(plane,1)==1
    line = [point n];
elseif size(point,1)==1
    line = [repmat(point, size(plane,1), 1) n];
elseif size(plane,1)==1
    line = [point repmat(n, size(point,1),1)];
else
    error('either 1-to-N, N-to-1, or N-to-N, but not N-to-M');
end
point = intersectPlaneLine(plane, line);