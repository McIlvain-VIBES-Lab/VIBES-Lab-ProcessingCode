function varargout = sph2cart2(varargin)
%SPH2CART2 cvonvert spherical coordinate to cartesian coordinate
%
%   usage :
%   C = SPH2CART2(S)
%   C = SPH2CART2(PHI, THETA)       (assume rho = 1)
%   C = SPH2CART2(PHI, THETA, RHO)   
%   [X, Y, Z] = SPH2CART2(PHI, THETA, RHO);
%
%   S = [phi theta rho] (sphercial coordiante).
%   C = [X Y Z]  (cartesian coordinate)
%
%   Math convention is used : theta is angle with vertical, 0 for north
%   pole, +pi for south pole, pi/2 for points with z=0.
%   phi is the same as matlab cart2sph : angle from Ox axis, counted
%   counter-clockwise.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY
%   22/03/2005 : make test for 2 args, and add radius if not specified for
%       1 arg.

if length(varargin)==1
    var = varargin{1};
    if size(var, 2)==2
        var = [var ones(size(var, 1), 1)];
    end
elseif length(varargin)==2
    var = [varargin{1} varargin{2} ones(size(varargin{1}))];
elseif length(varargin)==3
    var = [varargin{1} varargin{2} varargin{3}];
end

[x y z] = sph2cart(var(:,1), pi/2-var(:,2), var(:,3));

if nargout == 1 || nargout == 0
    varargout{1} = [x, y, z];
else
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
end
    