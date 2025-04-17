function varargout = drawPolyhedra(varargin)
%DRAWPOLYHEDRA : draw polyhedra defined by vertices and faces
%
%
%
%   H = DRAWPOLYHEDRA(...) also return handles to the created patches.
%
%
%   See also : drawPolygon
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 10/02/2005.
%


% process input parameters
if length(varargin)==2
    nodes = varargin{1};
    faces = varargin{2};
elseif length(varargin)==3
    nodes = varargin{1};
    faces = varargin{2};
    style = varargin{3};
else
    error ('wrong number of arguments in "drawPolyhedra"');
end

color = [1 0 0];

% --------------------
% main loop : for each face

hold on;
nf = size(faces, 1);
if iscell(faces)
    % array CELLS is a cell array
    for i=1:nf
        % get nodes of the cell
        cnodes = faces{i}';
        
        h(i) = patch(nodes(cnodes, 1), nodes(cnodes, 2), nodes(cnodes, 3), color);
    end
else
    
    % array FACES is a NC*NV indices array, with NV : number of vertices of
    % each face, and NC number of cells
    for i=1:nf
        % get nodes of the cell
        cnodes = faces(i,:)';
        
        h(i) = patch(nodes(cnodes, 1), nodes(cnodes, 2), nodes(cnodes, 3), color);
    end
end


% format output parameters
if nargout>0
    varargout{1}=h;
end