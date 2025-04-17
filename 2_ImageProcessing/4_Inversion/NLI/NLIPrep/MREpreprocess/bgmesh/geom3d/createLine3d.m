function line = createLine3d(varargin)
%CREATELINE3D create a line with various inputs.
%
%   Line is represented in a parametric form : [x0 y0 z0 dx dy dz]
%       x = x0 + t*dx
%       y = y0 + t*dy;
%       z = z0 + t*dz;
%
%
%   l = createLine3d(p1, p2) return the line going through the two given
%       points.
%   
%   l = createLine3d(x0, y0, z0, dx, dy, dz) the line going through point 
%       (x0, y0, z0) and with direction vector(dx, dy, dz).
%
%   l = createLine3d(P, dx, dy, dz) the line going through point P given by 
%       (x0, y0, z0) and with direction vector(dx dy dz).
%
%   
%   l = createLine3d(theta, phi) create a line originated at (0,0) and
%       with angles theta and phi.
%
%   l = createLine3d(rho, theta, phi) create a line with direction given by
%       theta and phi, and whose min distance to origin equals rho. rho can
%       be negative, in this case, the line is the same as with 
%       CREATELINE(-rho, theta+pi, phi+pi), but the orientation is
%       different. 
%
%
%   Note : in all cases, parameters can be vertical arrays of the same
%   dimension. The result is then an array of lines, of dimensions [N*4].
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%
%
% NOTE by Songbai JI:
%    The program actually does not support input argument of 3 or 4.  Need
%    to fix!

%   HISTORY :


%   NOTE : A 3d line can also be represented with a 1*7 array : 
%   [x0 y0 z0 dx dy dz t].
%   whith 't' being one of the following : 
%   - t=0 : line is a singleton (x0,y0)
%   - t=1 : line is an edge segment, between points (x0,y0) and (x0+dx,
%   y0+dy).
%   - t=Inf : line is a Ray, originated from (x0,y0) and going to infinity
%   in the direction(dx,dy).
%   - t=-Inf : line is a Ray, originated from (x0,y0) and going to infinity
%   in the direction(-dx,-dy).
%   - t=NaN : line is a real straight line, and contains all points
%   verifying the above equation.
%   This seems us a convenient way to represent uniformly all kind of lines
%   (including edges, rays, and even point).
%

if length(varargin)==1    
    error('Wrong number of arguments in ''createLine3d'' ');
    
elseif length(varargin)==2    
    % 2 input parameters. They can be :
    % - 2 points, then 2 arrays of 1*2 double.
    v1 = varargin{1};
    v2 = varargin{2};
    if size(v1, 2)==3
        % first input parameter is first point, and second input is the
        % second point.
        line = [v1(:,1) v1(:,2) v1(:,3) v2(:,1)-v1(:,1) v2(:,2)-v1(:,2) v2(:,3)-v1(:,3)];    
    end
    
elseif length(varargin)==3
    % 3 input parameters :
    error('Wrong number of arguments in ''createLine3d'' ');
   
elseif length(varargin)==4
    % 4 input parameters :
    error('Wrong number of arguments in ''createLine3d'' ');
else
    error('Wrong number of arguments in ''createLine3d'' ');
end
