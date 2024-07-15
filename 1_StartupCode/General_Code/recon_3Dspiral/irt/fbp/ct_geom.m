  function st = ct_geom(type, varargin)
%|function st = ct_geom(type, varargin)
%|
%| Create the "CT geometry" structure that describes the sampling
%| characteristics of a cone-beam CT system (axial or helical).
%| (Use sino_geom() for 2D fan-beam or parallel-beam systems.)
%|
%| in
%|	type	'fan' (multi-slice fan-beam)
%|
%| options for all geometries
%|	'orbit_start'		default: 0
%|	'orbit'			[degrees] default: 180 for parallel / mojette
%|					or 360 for fan (negative for CW)
%|	'down'			down-sampling factor, for testing
%|	'units'			string to print distance units (default: 'mm')
%|
%| options for fan-beam
%|	'ns'			# of horizontal samples
%|	'nt'			# of vertical samples
%|	'na' | 'nbeta'		# of angular samples
%|	'ds'			horizontal sample spacing (default: 1)
%|	'dt'			vertical sample spacing (default: -ds)
%|				or {'dz', dz} to use dz * dsd / dso (usual CT)
%|	'offset_s'		unitless fraction of a channel (default: 0)
%|				(relative to line between two center channels).
%|				use 0.25 or 1.25 for "quarter-detector offset"
%|	'offset_t'		unitless (default: 0)
%|
%| options for helical
%|	'pitch'			bed_travel_per_rotation / axial_fov. default: 0.
%|				(unitless. can be negative. usually near 1.0)
%|	'source_z0'		z-location of source for first view. Default 0.
%|					It must have same units as dt and ds.
%|
%|	'user_source_zs' [na]	user-specified source z-locations for each view.
%|				usually this is empty (default) in which case
%|				source_zs is computed internally from "pitch"
%|				It must have same units as dt and ds.
%|
%|	fan beam distances:
%|	'dsd' | 'dis_src_det'	default: inf (parallel beam)
%|	'dso' | 'dis_src_iso'	default: inf (parallel beam)
%|	'dod' | 'dis_iso_det'	default: 0
%|	'dfs' | 'dis_foc_src'	default: 0 (3rd generation CT arc),
%|					use 'inf' for flat detector
%|
%| out
%|	st	(struct)	initialized structure
%|
%| methods
%|	st.shape(sino)		reshape sinograms that are columns into 3d array
%|	st.s			s sample locations
%|	st.t			t sample locations
%|	st.ws			(ns-1)/2 + st.offset_s
%|	st.wt			(nt-1)/2 + st.offset_t
%|	st.ad			[na] source angles in degrees
%|	st.ar			[na] source angles in radians
%|	st.dim			dimensions: [st.ns st.nt st.na]
%|	st.downsample(down)	reduce sampling by integer factor
%|	st.ones			ones(ns,nt,na, 'single')
%|	st.zeros		zeros(ns,nt,na, 'single')
%|	st.rmax			max radius within FOV
%|	st.zfov			axial FOV
%|	st.source_zs		[na] z-locations of source for each view
%|	st.shape(sino(:))	reshape to [ns,nt,na,?]
%|	st.unitv(is,it,ia)	unit 'vector' with single nonzero element
%|	st.plot([ig])		show geometry
%|
%|	trick: you can make orbit=0 and orbit_start = column vector (length na)
%|	if you need nonuniformly spaced projection view angles.
%|
%| Copyright 2006-1-18, Jeff Fessler, University of Michigan
%|
%| 2009-12-04 modified source_zs definition to use source_z0, eliminate na/2

if nargin == 1 && streq(type, 'test'), ct_geom_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

if streq(type, 'ge1') % special case: GE fan-beam
	st = ct_geom_ge1(type, varargin{:});
return
end

if streq(type, 'ge2') % special case: GE axial or helical
	st = ct_geom_ge2(type, varargin{:});
return
end

% defaults
st.type = type;
st.ns = [];
st.nt = [];
st.na = [];
st.down = 1;
st.orbit_start = 0;
st.pitch = 0; % default for axial
st.source_z0 = 0; % default for axial
st.units = 'mm';
st.user_source_zs = [];

if streq(type, 'fan')
	st = ct_geom_fan(st, varargin{:});
%elseif streq(type, 'par')
%	st = ct_geom_par(st, varargin{:});
%elseif streq(type, 'moj')
%	st = ct_geom_moj(st, varargin{:});
else
	fail('unknown sinotype %s', type)
end

if isempty(st.na), st.na = 2 * floor(st.ns * pi/2 / 2); end

ct_geom_check_source_zs(st)

meth = { ...
	's', @ct_geom_s, '()'; ...
	't', @ct_geom_t, '()'; ...
	'ws', @ct_geom_ws, '()'; ...
	'wt', @ct_geom_wt, '()'; ...
	'ad', @ct_geom_ad, '()'; ...
	'ar', @ct_geom_ar, '()'; ...
	'zfov', @ct_geom_zfov, '()'; ...
	'source_dz_per_view', @ct_geom_source_dz_per_view, '()'; ...
	'source_zs', @ct_geom_source_zs, '() -> [na]'; ...
	'downsample', @ct_geom_downsample, '()'; ...
	'dim', @ct_geom_dim, '()'; ...
	'ones', @ct_geom_ones, '()'; ...
	'rmax', @ct_geom_rmax, '()'; ...
	'unitv', @ct_geom_unitv, '() | (is,it,ia)'; ...
	'zeros', @ct_geom_zeros, '()'; ...
	'shape', @ct_geom_shape, '()'; ...
	'plot', @ct_geom_plot, '() | (ig)';
	};

st = strum(st, meth);

if st.down ~= 1
	down = st.down; st.down = 1; % trick
	st = st.downsample(down);
end


% ct_geom_dim()
function dim = ct_geom_dim(st)
dim = [st.ns st.nt st.na];
if isempty(st.ns) || isempty(st.nt) || isempty(st.na)
	error 'dim requested without ns,nt,na'
end


% ct_geom_ones()
% sinogram of all ones
function out = ct_geom_ones(st)
out = ones(st.dim, 'single');


% ct_geom_zeros()
% sinogram of all zeros
function out = ct_geom_zeros(st)
out = zeros(st.dim, 'single');


% ct_geom_unitv()
% sinogram with a single ray
function out = ct_geom_unitv(st, is, it, ia)
out = st.zeros;
if ~isvar('is') || isempty(is)
	is = floor(st.ns/2 + 1);
	it = floor(st.nt/2 + 1);
	ia = 1;
end
out(is,it,ia) = 1;


%
% ct_geom_rmax()
%
% max radius within fov
function rmax = ct_geom_rmax(st)
smax = max(abs(st.s));
if streq(st.type, 'fan')
	if isinf(st.dso) % parallel
		rmax = smax;
	elseif st.dfs == 0 % arc
		rmax = st.dso * sin(smax / st.dsd);
	elseif isinf(st.dfs) % flat
		rmax = st.dso * sin(atan(smax / st.dsd));
	else
		error 'unknown case'
	end
end


% ct_geom_ws()
% 'middle' sample position
function ws = ct_geom_ws(st)
ws = (st.ns-1)/2 + st.offset_s;


% ct_geom_wt()
% 'middle' sample position
function wt = ct_geom_wt(st)
wt = (st.nt-1)/2 + st.offset_t;


% ct_geom_s()
% sample locations ('radial')
function s = ct_geom_s(st, varargin)
s = st.ds * ([0:st.ns-1]' - st.ws);
if length(varargin)
	s = s(varargin{:});
end


% ct_geom_t()
% sample locations ('radial')
function t = ct_geom_t(st, varargin)
t = st.dt * ([0:st.nt-1]' - st.wt);
if length(varargin)
	t = t(varargin{:});
end


% ct_geom_ad()
% angular sample locations (degrees)
function ang = ct_geom_ad(st, varargin)
ang = [0:st.na-1]'/st.na * st.orbit + st.orbit_start;
ang = ang(varargin{:});

% ct_geom_ar()
% angular sample locations (radians)
function ang = ct_geom_ar(st, varargin)
ang = deg2rad(ct_geom_ad(st));
ang = ang(varargin{:});


% ct_geom_shape()
% reshape into sinogram array
function sino = ct_geom_shape(st, sino)
sino = reshapee(sino, st.ns, st.nt, st.na, []);


%
% ct_geom_downsample()
%
% down-sample (for testing)
function st = ct_geom_downsample(st, down)
st.down = st.down * down;

if ~isempty(st.user_source_zs)
	st.user_source_zs = st.user_source_zs(1:down:st.na);
end

st.ns = 2 * round(st.ns / down / 2); % keep it even
st.nt = 2 * round(st.nt / down / 2); % keep it even
st.na = length([1:down:st.na]);

if streq(st.type, 'fan')
	st.ds = st.ds * down;
	st.dt = st.dt * down;
else
	error(['unknown sinotype ' type])
end


% ct_geom_zfov()
% axial FOV, considering magnification factor, as collimated at iso
function zfov = ct_geom_zfov(st)
zfov = st.dso / st.dsd * st.nt * st.dt; % collimated at iso


% ct_geom_source_dz_per_view()
% source axial travel for each view
function out = ct_geom_source_dz_per_view(st)
if ~isempty(st.user_source_zs), fail 'undefined', end
na_per_360 = st.na * (360 / st.orbit); % # views per turn
out = st.pitch * st.zfov / na_per_360;


% ct_geom_check_source_zs()
function ct_geom_check_source_zs(st)
if ~isempty(st.user_source_zs) && (st.pitch ~= 0)
	fail('only one of "pitch" (recommended) or user_source_zs can be given')
end
if ~isempty(st.user_source_zs) && length(st.user_source_zs) ~= st.na
	fail('user_source_zs size mismatch')
end


%
% ct_geom_source_zs()
%
% source z locations
function source_zs = ct_geom_source_zs(st, varargin)

if ~isempty(st.user_source_zs) % user-provided
	source_zs = st.user_source_zs;
	source_zs = source_zs(varargin{:});
else
	source_dz = st.source_dz_per_view;
	source_zs = st.source_z0 + [0:st.na-1]' * source_dz;
%	source_zs = (st.offset_z + [0:st.na-1]' - floor(st.na/2)) * source_dz;
	source_zs = source_zs(varargin{:});
end


%
% ct_geom_fan()
%
function st = ct_geom_fan(st, varargin);

% defaults
st.orbit = 360; % [degrees]
st.ds		= 1;
st.dt		= [];
st.offset_s	= 0;
st.offset_t	= 0;

st.dsd = [];	% dis_src_det
st.dso = [];	% dis_src_iso
st.dod = [];	% dis_iso_det
st.dfs = 0;	% dis_foc_src (3rd gen CT)

subs = { ...
	'src_det_dis', 'dsd';
	'dis_src_det', 'dsd';
	'dis_src_iso', 'dso';
	'dis_iso_det', 'dod';
	'dis_foc_src', 'dfs';
	'nbeta', 'na';
	'source_zs', 'user_source_zs';
	};
st = vararg_pair(st, varargin, 'subs', subs);

% work out distances
if (~isempty(st.dsd) && isinf(st.dsd)) ...
|| (~isempty(st.dso) && isinf(st.dso)) % handle parallel-beam case gracefully
	st.dsd = inf; st.dso = inf; st.dod = 1;
end
if isempty(st.dsd) + isempty(st.dso) + isempty(st.dod) > 1
	error 'must provide at least two of dsd, dso, dod'
end
if isempty(st.dsd), st.dsd = st.dso + st.dod; end
if isempty(st.dso), st.dso = st.dsd - st.dod; end
if isempty(st.dod), st.dod = st.dsd - st.dso; end
if st.dso + st.dod ~= st.dsd
	error 'bad fan-beam distances'
end

if isempty(st.dt), st.dt = -st.ds; end
if iscell(st.dt)
	if length(st.dt) == 2 && streq(st.dt{1}, 'dz')
		dz = st.dt{2};
		st.dt = dz * st.dsd / st.dso;
		printm('dt = dz * dsd / dso = %g * %g / %g = %g', ...
			dz, st.dsd, st.dso, st.dt)
	else
		fail 'bad dt cell'
	end
end


%
% ct_geom_plot2()
% picture of 2D source position / detector geometry
%
function ct_geom_plot2(st, ig)
if ~streq(st.type, 'fan'), error 'only fan done', end
x0 = 0;
y0 = st.dso;
t = linspace(0,2*pi,1001);

if isinf(st.dsd) % parallel beam
	rfov = max(abs(st.s));
	plot(	0, 0, '.', ...
		rfov * cos(t), rfov * sin(t), 'm:', ...
		st.s, -rfov, 'yo')
	

else % fan beam

	switch st.dfs
	case 0
		gam = st.s / st.dsd; % 3rd gen: equiangular
	case inf
		gam = atan(st.s / st.dsd); % flat
	otherwise
		error 'not done'
	end
	xds = st.dsd * sin(gam);
	yds = st.dso - st.dsd * cos(gam);
	rot = deg2rad(st.orbit_start);
	rot = [cos(rot) sin(rot); -sin(rot) cos(rot)];
	p0 = rot * [x0; y0];
	pd = rot * [xds'; yds'];
	rfov = st.dso * sin(max(abs(gam)));

	plot(	0, 0, '.', ...
		p0(1), p0(2), 's', ...
		[pd(1,1) p0(1) pd(1,end)], [pd(2,1) p0(2) pd(2,end)], '-', ...
		st.dso * cos(t), st.dso * sin(t), '--', ... % source circle
		rfov * cos(t), rfov * sin(t), 'm:', ...
		pd(1,:), pd(2,:), 'yo')
end


if isvar('ig') && ~isempty(ig)
	hold on
	xmin = min(ig.x); xmax = max(ig.x);
	ymin = min(ig.y); ymax = max(ig.y);
	plot([xmax xmin xmin xmax xmax], [ymax ymax ymin ymin ymax], 'g-')
	hold off
end
title(sprintf('fov = %g', rfov))
axis square, zoom on


%
% ct_geom_plot3()
% picture of 3D helical CT geometry
%
function ct_geom_plot3(st, ig)
if ~streq(st.type, 'fan'), error 'only fan done', end

t1 = -st.dso * sin(st.ar); % src x pos
t2 = st.dso * cos(st.ar); % src r pos
t3 = st.source_zs; % src z pos
plot3(t1, t2, t3, 'y.') % trajectory
view(-110, 22)
hold on
text(t1(1), t2(1), t3(1), '1')
center = 1 + floor(st.na/2);
text(t1(center), t2(center), t3(center), sprintf('%d',center))
text(t1(end), t2(end), t3(end), sprintf('%d',st.na))

% axes
plot3([-100*ceil(max(abs(t1))/100) 100*ceil(max(abs(t1))/100)], [0 0], [0 0])
text(100*ceil(max(abs(t1))/100),0, 'x')
plot3([0 0], [-100*ceil(max(abs(t2))/100) 100*ceil(max(abs(t2))/100)], [0 0])
text(0, 100*ceil(max(abs(t2))/100), 'y')
plot3([0 0], [0 0], [-10*ceil(max(abs(t3))/10) 10*ceil(max(abs(t3))/10)])
text(0, 0, 10*ceil(max(abs(t3))/10), 'z')

grid on
xlabelf('x (%s)', st.units)
ylabelf('y (%s)', st.units)
zlabelf('z (%s)', st.units)

if isvar('ig') && ~isempty(ig) && isfield(ig.meth,'z')
	hold on
	xmin = min(ig.x); xmax = max(ig.x);
	ymin = min(ig.y); ymax = max(ig.y);
	zmin = min(ig.z); zmax = max(ig.z);
	plot3(	[xmin xmax xmax xmin xmin ], ...
		[ymin ymin ymax ymax ymin ], ...
		[zmin zmin zmin zmin zmin ], ...
		'g-')
	plot3(	[xmin xmax xmax xmin xmin], ...
		[ymin ymin ymax ymax ymin], ...
		[zmax zmax zmax zmax zmax], ...
		'g-')
	plot3(xmin*[1 1], ymin*[1 1], [zmin zmax], 'g-')
	plot3(xmax*[1 1], ymin*[1 1], [zmin zmax], 'g-')
	plot3(xmin*[1 1], ymax*[1 1], [zmin zmax], 'g-')
	plot3(xmax*[1 1], ymax*[1 1], [zmin zmax], 'g-')
end

% detector array for source at first position
cen_loc = [t1(1) t2(1) t3(1)];
unit_vec = [-cen_loc(1) -cen_loc(2) 0] / sqrt(cen_loc(1)^2+cen_loc(2)^2);
% diametrically opposite point of the source on the detector
det_cen_loc = cen_loc + st.dsd * unit_vec;
plot3(det_cen_loc(1),det_cen_loc(2),det_cen_loc(3),'rs')
theta = st.s/st.dsd;
det_x = cen_loc(1) + st.dsd * (unit_vec(1)*cos(theta)-unit_vec(2)*sin(theta));
det_y = cen_loc(2) + st.dsd * (unit_vec(1)*sin(theta)+unit_vec(2)* ...
				cos(theta));
st_t = st.t;
for tt=1:length(st_t)
	plot3(det_x, det_y, cen_loc(3)+st_t(tt)*ones(size(theta)), 'r.')
end

% detector array for source at middle position
cen_loc = [t1(center) t2(center) t3(center)];
unit_vec = [-cen_loc(1) -cen_loc(2) 0] / sqrt(cen_loc(1)^2+cen_loc(2)^2);

% diametrically opposite point of the source on the detector
det_cen_loc = cen_loc + st.dsd * unit_vec;
plot3(det_cen_loc(1),det_cen_loc(2),det_cen_loc(3),'rs')
theta = st.s/st.dsd;
det_x = cen_loc(1) + st.dsd * (unit_vec(1)*cos(theta)-unit_vec(2)*sin(theta));
det_y = cen_loc(2) + st.dsd * (unit_vec(1)*sin(theta)+unit_vec(2)* ...
				cos(theta));
st_t = st.t;
for tt=1:length(st_t)
	plot3(det_x, det_y, cen_loc(3)+st_t(tt)*ones(size(theta)), 'c.')
end

% detector array for source at final position
cen_loc = [t1(end) t2(end) t3(end)];
unit_vec = [-cen_loc(1) -cen_loc(2) 0] / sqrt(cen_loc(1)^2+cen_loc(2)^2);
% diametrically opposite point of the source on the detector
det_cen_loc = cen_loc + st.dsd * unit_vec;
plot3(det_cen_loc(1),det_cen_loc(2),det_cen_loc(3),'rs')
theta = st.s/st.dsd;
det_x = cen_loc(1) + st.dsd * (unit_vec(1)*cos(theta)-unit_vec(2)*sin(theta));
det_y = cen_loc(2) + st.dsd * (unit_vec(1)*sin(theta)+unit_vec(2)* ...
				cos(theta));
st_t = st.t;
for tt=1:length(st_t)
	plot3(det_x, det_y, cen_loc(3)+st_t(tt)*ones(size(theta)), 'm.')
end


%
% ct_geom_plot()
%
function out = ct_geom_plot(st, varargin)
im clf
if all(st.source_zs == 0) % 2D or axial case
	ct_geom_plot2(st, varargin{:});
else
	ct_geom_plot3(st, varargin{:});
end
if nargout, out = []; end


%
% ct_geom_ge1()
% 'lightspeed';
% these numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl
%
function st = ct_geom_ge1(type, varargin)
st = ct_geom('fan', ...
	'ns', 888, ... % detector channels
	'nt', 1, ... % detector rows
	'na', 984, ... % angular samples
	'orbit', 360, ...
	'offset_s', 1.25, ... % quarter-detector offset
	'dsd', 949.075, ...
	'dod', 408.075, ...
	'dfs', 0, ... % arc
	'ds', 1.0239, ... % detector pitch
	'dt', 1.0964, ... % detector row spacing for 0.625mm slices, 2009-12-06
	varargin{:});
% 'strip_width', [], ...

%
% ct_geom_ge2()
% helical CT
%
function st = ct_geom_ge2(type, varargin)
st = ct_geom_ge1(type, ...
	'nt', 64, ... % 64-slice
	'dt', 949.0750 / 541 * 0.625, ... % about 1.0964
	varargin{:});


%
% ct_geom_test()
%
function ct_geom_test
% axial cone-beam test
cg = ct_geom('fan', 'ns', 888, 'nt', 64, 'na', 984, ...
	'offset_s', 1.25, ...
	'dsd', 949, 'dod', 408);
cg.ad(2);
cg.downsample(2);
cg.rmax;
cg.ws;
cg.s(cg.ns/2+1);
if im
	cg.plot;
prompt
end

% test user_source_zs
cg = ct_geom('fan', 'dsd', 949, 'dod', 408, ...
	'ns', 888, 'ds', 1.0239, 'offset_s', 1.25, ...
	'nt', 8, 'dt', 1.0964, ...
	'na', 2*984, 'orbit', 2*360, ...
	'source_zs', []);
cg.source_zs(1:2:7);

% helical cone-beam test
cg = ct_geom('fan', 'dsd', 949, 'dod', 408, ...
	'ns', 888, 'ds', 1.0239, 'offset_s', 1.25, ...
	'nt', 8, 'dt', 1.0964, ...
	'na', 2*984, 'orbit', 2*360, ...
	'pitch', 1.0);

cg.source_zs(1:2:7);
cg.downsample(2);
cg.s(cg.ns/2+1);
if im
	cg.plot;
prompt
end
