function alpha = sphericalAngle(p1, p2, p3)
%SPHERICALANGLE compute angle on the sphere
%
%   ALPHA = SPHERICALANGLE(P1, P2, P3)
%   compute angle (P1, P2, P2), in radians, between 0 and 2*PI.
%
%   Points are given either as [x y z] (there will be normalized to lie on
%   the unit sphere), or as [phi theta], with phi being the longitude in [0
%   2*PI] and theta being the elevation on horizontal [-pi/2 pi/2].
%
%
%   TODO :
%   - lame implementation, give result only <PI.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY

% test if points are given as matlab spherical coordinate
if size(p1, 2) ==2
    [x y z] = sph2cart(p1(:,1), p1(:,2));
    p1 = [x y z];
    [x y z] = sph2cart(p2(:,1), p2(:,2));
    p2 = [x y z];
    [x y z] = sph2cart(p3(:,1), p3(:,2));
    p3 = [x y z];
end

% normalize points
p1 = normalize(p1);
p2 = normalize(p2);
p3 = normalize(p3);


a12 = anglePoints3d(p1, [0 0 0], p2);
a23 = anglePoints3d(p2, [0 0 0], p3);

alpha = pi + real(asin(det([p2;p3;p1])/sin(a12)/sin(a23)));

plane = createPlane(p2, p2);

pi1 = planePosition(intersectPlaneLine(plane, [0 0 0 p1]), plane);
pi3 = planePosition(intersectPlaneLine(plane, [0 0 0 p3]), plane);
alpha = angle3Points(pi1, [0 0], pi3);