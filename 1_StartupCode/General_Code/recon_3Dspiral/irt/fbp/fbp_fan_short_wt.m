  function wt = fbp_fan_short_wt(sg, varargin)
%|function wt = fbp_fan_short_wt(sg, [options])
%| Sinogram weighting for fan-beam short scan
%| in
%|	sg	strum		sino_geom
%| option
%|	type	'parker'	from parker:82:oss (Med Phys 1982)
%| out
%|	wt	[nb na]
%|
%| Copyright 2009-12-10, Jeff Fessler and Janghwan Cho, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(sg, 'test'), fbp_fan_short_wt_test, return, end

arg.type = 'parker';
arg = vararg_pair(arg, varargin);

switch arg.type
case 'parker'
	wt = fbp_fan_short_wt_parker(sg);
otherwise
	fail('unknown type %s', arg.type)
end


% fbp_fan_short_wt_parker()
function wt = fbp_fan_short_wt_parker(sg)
nb = sg.nb;
na = sg.na;
bet = sg.ar;
gam = sg.gamma;
[g b] = ndgrid(gam, bet);
gammax = sg.gamma_max; % half of fan angle

fun = @(x) sin(pi/2 * x).^2; % smooth out [0,1] ramp
% todo: could use integrals of this function over the
% tiny angular range of each projection view so that
% the sum over beta of these functions is a constant.

%wt = nan(nb,na);
wt = ones(nb,na);
ii = 0 <= b & b <  2 * (gammax - g);
tmp = b(ii) ./ (2 * (gammax - g(ii)));
wt(ii) = fun(tmp);

%ii = 2 * (gammax - g) < b & b < pi - 2 * g;
%wt(ii) = 1;

ii = pi - 2 * g < b & b <= pi + 2 * gammax;
tmp = (pi + 2*gammax - b(ii)) ./ (2 * (gammax + g(ii)));
wt(ii) = fun(tmp);

%if any(isnan(wt(:))), warn 'nan!', end


% fbp_fan_short_wt_test
function fbp_fan_short_wt_test

% generate object, image geometry
ig = image_geom('nx', 128, 'dx', 4);
ell = [	0 0 50 40 0 1;
	-18 0 12 12 0 1;
	18 -0 12 12 0 1] * 4;
xtrue = ellipse_im(ig, ell, 'oversample', 4);

% sinogram geometry for both short and full scan
sg_short = sino_geom('ge1', 'orbit', 'short', 'down', 4);
sg_360 = sino_geom('ge1', 'down', 4);

% sinogram of the object for each scan
sino_short = ellipse_sino(sg_short, ell, 'oversample', 2);
sino_360 = ellipse_sino(sg_360, ell, 'oversample', 2);

% apply parker weighting
wt = fbp_fan_short_wt(sg_short); % parker weighting
scale = sg_short.orbit / 180; % scale factor due to orbit
sino_parker = sino_short .* wt * scale;

% FBP reconstructed images
fbp_geom_short = fbp2(sg_short, ig);
fbp_geom_360 = fbp2(sg_360, ig);

fbp_w_short = fbp2(sino_short, fbp_geom_short);
fbp_w_parker = fbp2(sino_parker, fbp_geom_short);
fbp_w_360 = fbp2(sino_360, fbp_geom_360);

% plot
im plc 4 3
im subplot [1 2 3]
im(rad2deg(sg_short.gamma), sg_short.ad, wt), cbar
xlabel 'gamma [degrees]'
ylabel 'beta [degrees]'
title 'Parker weighting'
clim = [0 8];
im(4, fbp_w_360, clim, 'FBP: full scan'), cbar
im(5, fbp_w_short, clim, 'FBP: short scan w/o parker weighting'), cbar
im(6, fbp_w_parker, clim, 'FBP, short scan w/ parker weighting'), cbar
im(7, xtrue, clim, 'True'), cbar
im subplot [8 9]
plot([fbp_w_360(:,end/2) fbp_w_short(:,end/2) fbp_w_parker(:,end/2)])
title 'middle slice, y = end/2'
legend('full', 'short w/o parker', 'short w/ parker')

im(10, sino_360)
im(11, sino_short)
im(12, sino_parker)

%yaxis_pi '0 p'
%plot(diff(wt,1))
%savefig fig_tomo_fan_short_wt
