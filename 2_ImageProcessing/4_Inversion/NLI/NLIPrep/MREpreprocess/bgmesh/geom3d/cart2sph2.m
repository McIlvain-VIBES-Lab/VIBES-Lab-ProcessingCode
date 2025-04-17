function varargout = cart2sph2(varargin)
%CART2SPH2 convert cartesian 2 spherical coordinate
%
%   usage :
%   S = CART2SPH2(C)
%   C = [X Y Z]  (cartesian coordinate
%   S = [phi theta rho] (sphercial coordiante).
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

if length(varargin)==1
    var = varargin{1};
elseif length(varargin)==3
    var = [varargin{1} varargin{2} varargin{3}];
end

[t p r] = cart2sph(var(:,1), var(:,2), var(:,3));

if nargout == 1 || nargout == 0
    varargout{1} = [t pi/2-p r];
else
    varargout{1} = t;
    varargout{2} = pi/2-p;
    varargout{3} = r;
end
    