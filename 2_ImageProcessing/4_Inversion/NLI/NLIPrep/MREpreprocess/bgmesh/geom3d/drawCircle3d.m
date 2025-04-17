function varargout = drawCircle3d(varargin)
%DRAWCIRCLE3D draw a 3D circle
%
%   Ps=ossible calls for the function :
%   DRAWCIRCLE3D([XC YC ZC R PHI THETA])            1
%   DRAWCIRCLE3D([XC YC ZC R PHI THETA PSI])        1
%   DRAWCIRCLE3D([XC YC ZC R], [PHI THETA])         2
%   DRAWCIRCLE3D([XC YC ZC R], [PHI THETA PSI])     2
%   DRAWCIRCLE3D([XC YC ZC R], PHI, THETA)          3
%   DRAWCIRCLE3D([XC YC ZC], R, PHI, THETA)         4
%   DRAWCIRCLE3D([XC YC ZC R], PHI, THETA, PSI)     4
%   DRAWCIRCLE3D([XC YC ZC], R, PHI, THETA, PSI)    5
%   DRAWCIRCLE3D(XC, YC, ZC, R, PHI, THETA)         6
%   DRAWCIRCLE3D(XC, YC, ZC, R, PHI, THETA, PSI)    7
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005
%

%   HISTORY

if length(varargin)==1
    % get center and radius
    circle = varargin{1};
    xc = circle(:,1);
    yc = circle(:,2);
    zc = circle(:,3);
    r  = circle(:,4);
    
    % get angle of normal
    phi     = circle(:,5);
    theta   = circle(:,6);

    % get roll
    if size(circle, 2)==7
        psi = circle(:,7);
    else
        psi = zeros(size(circle, 1), 1);
    end
    
elseif length(varargin)==2
    % get center and radius
    sphere = varargin{1};
    xc = sphere(:,1);
    yc = sphere(:,2);
    zc = sphere(:,3);
    r  = sphere(:,4);
    
    % get angle of normal
    angle = varargin{2};
    phi     = angle(:,1);
    theta   = angle(:,2);
    
    % get roll
    if size(angle, 2)==3
        psi = angle(:,3);
    else
        psi = zeros(size(angle, 1), 1);
    end

elseif length(varargin)==3
    
    % get center and radius
    sphere = varargin{1};
    xc = sphere(:,1);
    yc = sphere(:,2);
    zc = sphere(:,3);
    r  = sphere(:,4);
    
    % get angle of normal and roll
    phi     = varargin{2};
    theta   = varargin{3};
    psi     = zeros(size(phi, 1), 1);
    
elseif length(varargin)==4
    % get center and radius
    sphere = varargin{1};
    xc = sphere(:,1);
    yc = sphere(:,2);
    zc = sphere(:,3);
    
    if size(sphere, 2)==4
        r   = sphere(:,4);
        phi     = varargin{2};
        theta   = varargin{3};
        psi     = varargin{4};
    else
        r   = varargin{2};
        phi     = varargin{3};
        theta   = varargin{4};
        psi     = zeros(size(phi, 1), 1);
    end
    
elseif length(varargin)==5
    % get center and radius
    sphere = varargin{1};
    xc = sphere(:,1);
    yc = sphere(:,2);
    zc = sphere(:,3);
    r  = varargin{2};
    phi     = varargin{3};
    theta   = varargin{4};
    psi     = varargin{5};

elseif length(varargin)==6
    xc      = varargin{1};
    yc      = varargin{2};
    zc      = varargin{3};
    r       = varargin{4};
    phi     = varargin{5};
    theta   = varargin{6};
    psi     = zeros(size(phi, 1), 1);
  
elseif length(varargin)==7
   
    xc      = varargin{1};
    yc      = varargin{2};
    zc      = varargin{3};
    r       = varargin{4};
    phi     = varargin{5};
    theta   = varargin{6};
    psi     = varargin{7};

else
    error('DRAWCIRCLE3D : please specify center and radius');
end

N = 64;
t = [0:2*pi/N:2*pi*(1-1/N) 2*pi];


x = r*cos(t)';
y = r*sin(t)';
z = zeros(length(t), 1);

circle0 = [x y z];

tr = translation3d(xc, yc, zc);
rot1 = rotationOz(psi);
rot2 = rotationOy(-theta);
rot3 = rotationOz(-phi);
trans = tr*rot3*rot2*rot1;

circle = transformPoint3d(circle0, trans);

h = drawCurve3d(circle);


if nargout>0
    varargout{1}=h;
end

