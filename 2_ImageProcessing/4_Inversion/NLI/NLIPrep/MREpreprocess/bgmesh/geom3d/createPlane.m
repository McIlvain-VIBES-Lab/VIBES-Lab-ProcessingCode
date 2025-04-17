function plane = createPlane(varargin)
%CREATEPPLANE create a plane in parametrized form
%
%   Create a plane in the following format : 
%   [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], where :
%   - (X0, Y0, Z0) is a point belonging to the plane
%   - (DX1, DY1, DZ1) is a first direction vector
%   - (DX2, DY2, DZ2) is a second direction vector
%   
%
%
%   PL = CREATEPLANE(P1, P2, P3) create a plane containing the 3 points
%
%   PL = CREATEPLANE(P0, N) create a plane from a point and from a normal
%   to the plane.
%   
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY :

if length(varargin)==2
    
    p0 = varargin{1};
    
    var = varargin{2};
    if size(var, 2)==2
        n = sph2cart2([var repmat(1, [size(var, 1) 1])]);
    elseif size(var, 2)==3
        n  = normalize(var);
    else
        error ('wrong number of parameters in createPlane');
    end
    
    % find a vector not colinear to the normal
    v0 = repmat([1 0 0], [size(p0, 1) 1]);    
    if abs(cross(n, v0, 2))<1e-14
        v0 = repmat([0 1 0], [size(p0, 1) 1]);
    end
    
    % create direction vectors
    v1 = normalize(cross(n, v0, 2));
    v2 = -normalize(cross(v1, n, 2));
    
    plane = [p0 v1 v2];
    return;
    
elseif length(varargin)==3
    p1 = varargin{1};    
    p2 = varargin{2};
    p3 = varargin{3};
    
    % create direction vectors
    v1 = p2-p1;
    v2 = p3-p1;
   
    plane = [p1 v1 v2];
    return;
  
else
    error('wrong number of arguments in "createPlane".');
end
    



    
    
    