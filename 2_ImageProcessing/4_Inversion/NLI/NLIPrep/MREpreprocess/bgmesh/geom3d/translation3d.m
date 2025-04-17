function trans = translation(varargin)
%TRANSLATION3D : return 4x4 matrix of a 3D translation
%
%   usage :
%   TRANS = TRANSLATION(DX, DY);
%   return the translation corresponding to DX and DY.
%   The returned matrix has the form :
%   [1 0 DX]
%   [0 1 DY]
%   [0 0  1]
%
%   TRANS = TRANSLATION(POINT);
%   return the translation corresponding to the given point [x y].
%
%
%   See also :
%   transformPoint, rotation
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 06/04/2004.
%

if length(varargin)==0
    dx = 0;
    dy = 0;
elseif length(varargin)==1
    var = varargin{1};
    dx = var(1);
    dy = var(2);
    dz = var(3);
else
    dx = varargin{1};
    dy = varargin{2};
    dz = varargin{3};
end

trans = [1 0 0 dx ; 0 1 0 dy ; 0 0 1 dz; 0 0 0 1];