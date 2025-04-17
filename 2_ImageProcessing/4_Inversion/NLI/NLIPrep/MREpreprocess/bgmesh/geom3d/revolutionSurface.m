function varargout = revolutionSurface(varargin)
%REVOLUTIONSURFACE create a surface of revolution from a planar curve
%
%   usage : 
%   [X Y Z] = revolutionSurface(XT, YT, N);
%   create the surface of revolution of parametrized function (xt, yt),
%   with N equally spaced slices.
%
%   [X Y Z] = revolutionSurface(CURVE, N);
%   is the same, but generating curve is given in a single parameter CURVE,
%   which is a [Nx2] array of 2D points.
%
%   Surface can be displayed using :
%   H = surf(X, Y, Z);
%   H is a handle to the created patch.
%
%   revolutionSurface(...);
%   by itself, directly shows the created patch.
%
%
%   TODO : add possibility to specify axis of revolution
%
%
%   ------
%   Author: David Legland
%   e-mail: david.legland@jouy.inra.fr
%   Created: 2004-04-09
%   Copyright 2005 INRA - CEPIA Nantes - MIAJ Jouy-en-Josas.

%   based on function cylinder from matlab

if length(varargin)==2
    curve = varargin{1};
    xt = curve(:,1);
    yt = curve(:,2);
    n = varargin{2};
elseif length(varargin)==3
    xt = varargin{1};
    xt = varargin{2};
    n = varargin{3};
end

% ensure length is enough
m = length(xt);
if m==1
    xt = [xt xt];
    m = 2;
end
        
% create revolution angles
theta = (0:n)/n*2*pi;
sintheta = sin(theta); sintheta(n+1) = 0;

% compute surface vertices
x = yt * cos(theta);
y = yt * sintheta;
z = xt * ones(size(theta));

% format output depending on how many output parameters
if nargout == 0
    surf(x,y,z)
elseif nargout==3
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
end


