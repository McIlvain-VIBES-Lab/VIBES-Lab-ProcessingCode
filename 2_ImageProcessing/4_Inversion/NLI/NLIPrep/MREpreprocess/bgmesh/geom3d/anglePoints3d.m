function alpha = anglePoints3d(varargin)
%ANGLEPOINTS3D compute angle between 2 3D points
%
%   ALPHA = ANGLEPOINTS3D(P1, P2)
%   compute angle (P1, O, P2), in radians, between 0 and PI.
%
%   ALPHA = ANGLEPOINTS3D(P1, P2, P3)
%   compute angle (P1, P2, P3), in radians, between 0 and PI.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY


p2 = [0 0 0];
if length(varargin)==2
    p1 = varargin{1};
    p0 = [0 0 0];;
    p2 = varargin{2};
elseif length(varargin)==3
    p1 = varargin{1};
    p0 = varargin{2};
    p2 = varargin{3};
end


% normalized points
p1 = normalize(p1-p0);
p2 = normalize(p2-p0);
alpha = acos(dot(p1, p2, 2));
