% ct_fan_beam_example.m
% compare FBP and iterative reconstruction for a 2D fan-beam CT problem
% Copyright 2009-10-05, Jeff Fessler, University of Michigan

if ~isvar('A'), printm 'setup geometry, image, sinogram'
	down = 4;
	ig = image_geom('nx', 512, 'fov', 50, 'down', down);
	sg = sino_geom('ge1', 'units', 'cm', 'down', down);

	% read image
	ddir = path_find_dir([filesep 'data']);
	xtrue256 = fld_read([ddir filesep 'ncat,256,slice,140,ct,x100.fld']);
	xtrue256 = xtrue256 / 200 * 0.4; % convert to 1/cm units 

	if 1 % more realistic sinogram from finer image
		ig_big = image_geom('nx', 512, 'fov', ig.fov, 'down', 2);
		Abig = Gtomo2_dscmex(sg, ig_big);
		sino_true = Abig * xtrue256;
	end
	xtrue = downsample2(xtrue256, 2);

	% system object
	A = Gtomo2_dscmex(sg, ig);

	im clf, im pl 2 2
	clim = [0 0.4];
	im(1, xtrue, 'x', clim), cbar
	im(2, sino_true, 'sino'), cbar

	clear ddir ig_big Abig
prompt
end


if ~isvar('sino'), printm 'noisy fan-beam data'
	I0 = 1e5; % incident photons
	rand('state', 0)
	% transmission data:
	yi = poisson(I0 * exp(-sino_true), 0, 'factor', 0.4); % poisson noise
	if any(yi(:) == 0)
		warn('%d of %d values are 0 in sinogram!', ...
			sum(yi(:)==0), length(yi(:)));
	end
	sino = log(I0 ./ max(yi,1)); % noisy fan-beam sinogram
	im(4, sino, 'noisy sino'), cbar
prompt
end


if ~isvar('fbp'), printm 'fbp 2d fan-beam reconstruction'
	tmp = fbp2(sg, ig);
	fbp = fbp2(sino, tmp);
	im(3, fbp, 'FBP', clim), cbar
return
end


% todo: iterative reconstruction
