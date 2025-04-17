function varargout = drawCylinder(cyl, varargin)
%DRAWCYLINDER draw a cylinder
%
%   usage :
%   drawCylinder(CYL)
%   where CYL is a cylinder defined by [x1 y2 z1 x2 y2 z2 r], from starting
%   and ending point of cylinder, together with radius, draws the
%   corresponding shape on the current axis.
%
%   drawCylinder(CYL, N)
%   uses N points for discretisation of angle
%
%   
%   drawCylinder(..., OPT)
%   with OPT = 'open' or 'closed', specify if bases of cylinder should be
%   drawn.
%
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/09/2005
%

%   HISTORY

if iscell(cyl)
    for i=1:length(cyl)
        res(i) = drawCylinder(cyl{i}, varargin{:});
    end
    
    if nargout>0
        varargout{1} = res;
    end
    
    return;
end

% default
N = 32;
closed = false;


% process input arguments
if length(varargin)==1
    var = varargin{1};
    if ischar(var)
        if strcmp(var, 'open')
            closed = false;
        elseif strcmp(var, 'closed')
            closed = true;
        end
    else
        N = var;
    end
elseif length(varargin)==2
    N = varargin{1};
    var = varargin{2};
   if strcmp(var, 'open')
        closed = false;
    elseif strcmp(var, 'closed')
        closed = true;
   end

end



[phi theta rho] = cart2sph2(cyl(4:6)-cyl(1:3));
dphi = 0:2*pi/N:2*pi;

r = cyl(7);

x = repmat(cos(dphi)*r, [2 1]);
y = repmat(sin(dphi)*r, [2 1]);
z = repmat([0;rho], [1 length(dphi)]);

pts = transformPoint3d([x(:) y(:) z(:)], rotationOy(-theta));
pts = transformPoint3d(pts, rotationOz(-phi));
pts = transformPoint3d(pts, translation3d(cyl(1:3)));

x2 = reshape(pts(:,1), size(x));
y2 = reshape(pts(:,2), size(x));
z2 = reshape(pts(:,3), size(x));


if nargout == 0
    surf(x2,y2,z2, 'FaceColor', 'g', 'edgeColor', 'none');
else
    xx = x2; yy = y2; zz = z2;
end
