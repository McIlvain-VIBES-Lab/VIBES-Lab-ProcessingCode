function varargout = drawCurve3d(varargin)
%DRAWCURVE3D Draw a 3D curve specified by a list of points
%
%   DRAWCURVE3D(COORD) packs coordinates in a single [N*3] array.
%
%   DRAWCURVE(PX, PY, PZ) specify coordinates in separate arrays.
%
%   H = DRAWCURVE(...) also return a handle to the list of line objects.
%
%   See Also :
%   drawPolygon
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%



px = 0;
py = 0;
pz = 0;
closed = false;
symbol = 'b-';

% check case we xant to draw several curves, stroed in a cell array
if length(varargin)>0
    var = varargin{1};
    if iscell(var)
        hold on;
        for i=1:length(var)
            h(i) = drawCurve3d(var{i}, varargin{2:end});
        end
    end
end


if length(varargin)==1
    var = varargin{1};
    if iscell(var)
        for i=1:length(var)
            h(i) = drawCurve3d(var{i});
        end
    else
        px = var(:, 1);
        py = var(:, 2);
        pz = var(:, 3);
    end
elseif length(varargin)==2
    var = varargin{1};
    px = var(:,1);
    py = var(:,2);
    pz = var(:,3);
    
    if ischar(varargin{2})
        
        if strmatch('close', varargin{2})
            closed = true;
        elseif strmatch('open', varargin{2})
            closed = false;
        else
            symbol = varargin{2};
        end               
    end
elseif length(varargin)==3
    px = varargin{1};
    py = varargin{2};
    pz = varargin{3};

elseif length(varargin)>=4
    px = varargin{1};
    py = varargin{2};
    pz = varargin{3};
    if strmatch('close', varargin{3})
        closed = true;
    elseif strcmp('open', varargin{3})
        closed = false;
    else
        symbol = varargin{3};
    end               
else
    error ('wrong number of arguments in "drawCurve3d"');
end

if closed
    px = [px; px(1)];
    py = [py; py(1)];
    pz = [pz; pz(1)];
end
h = plot3(px, py, pz, symbol);

if nargout>0
    varargout{1}=h;
end