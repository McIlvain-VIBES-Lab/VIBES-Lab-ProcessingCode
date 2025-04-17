function varargout = drawPoint3d(varargin)
%DRAWPOINT3D draw 3D point on the current axis.
%
%   DRAWPOINT(X, Y, Z) will draw points defined by coordinates X and Y.
%   X and Y are N*1 array, with N being number of points to be drawn.
%   If coordinates of points lie outside the visible area, points are
%   not drawn.
%
%   DRAWPOINT(COORD) packs coordinates in a single [N*3] array.
%
%   DRAWPOINT(..., OPT) will draw each point with given option. OPT is a
%   string compatible with 'plot' model. OPT is a single string, it is not
%   a string array.
%
%
%   H = DRAWPOINT(...) also return a handle to each of the drawn points.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY


px = 0;
py = 0;
pz = 0;
symbol = 'o';

if length(varargin)==1
    var = varargin{1};
    px = var(:, 1);
    py = var(:, 2);
    pz = var(:, 3);
elseif length(varargin)==2
    var = varargin{1};
    px = var(:,1);
    py = var(:,2);
    pz = var(:,3);
    symbol = varargin{2};
elseif length(varargin)==3
    px = varargin{1};
    py = varargin{2};
    pz = varargin{3};
elseif length(varargin)==4
    px = varargin{1};
    py = varargin{2};
    pz = varargin{3};
    symbol = varargin{4};
else
    error ('wrong number of arguments in "drawPoint"');
end


lim = get(gca, 'xlim');
xmin = lim(1);
xmax = lim(2);
lim = get(gca, 'ylim');
ymin = lim(1);
ymax = lim(2);
lim = get(gca, 'zlim');
zmin = lim(1);
zmax = lim(2);

% check validity for display
ok = px>=xmin;
ok = ok & px<=xmax;
ok = ok & py>=ymin;
ok = ok & py<=ymax;
ok = ok & pz>=zmin;
ok = ok & pz<=zmax;

h = plot3(px(ok), py(ok), pz(ok), symbol);

if nargout>0
    varargout{1}=h;
end