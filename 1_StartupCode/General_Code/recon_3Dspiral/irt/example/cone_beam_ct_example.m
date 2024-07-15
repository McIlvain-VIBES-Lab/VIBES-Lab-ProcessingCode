% cone_beam_ct_example.m
% Illustrate iterative cone-beam X-ray CT image reconstruction.
% This illustration uses a tiny system size because my cone-beam
% projector is very slow.  Users with their own fast cone-beam
% projector/backprojector can use this as a guide.
% Hopefully someday I will have faster cone-beam code...
%
% Copyright 2005-1-21, Jeff Fessler, University of Michigan

% First run the FDK example.  It generates the true image xtrue
% and noiseless projection views "proj" and noisy data "yi"
% and generates (noisy) FDK recon "xfdk" for comparison / initialization.
if ~isvar('xfdk')
	bi = 1e6; % 1M photons / ray
	ri = 0; % no scatter etc. for now
	dfs = inf; % flat!
	feldkamp_example
prompt
end

if 0, printm 'dd check' % check distance-driven vs analytical
	Ad = Gtomo_dd(cg, ig);
	pd = Ad * xtrue;
	im clf, im_toggle(proj(:,:,1:12:end), pd(:,:,1:12:end), [0 4.4])
	nrms(pd, proj)
return
end

% 3l system matrix
if ~isvar('A'), printm 'A'
	if 1
		A = Gcone(cg, ig, 'type', 'sf1');
	elseif exist('dd_ge2_mex') == 3 % for UM only!
		A = Gcone(cg, ig, 'type', 'dd2');
	else
		f.sys_type = aspire_pair(cg, ig, 'system', '3l');
		A = Gtomo3(f.sys_type, ig.mask, ig.nx, ig.ny, ig.nz, ...
			'chat', 0, 'permute213', true, 'checkmask', im&0);
	end
end

% block object for ordered-subsets iterations
if ~isvar('Ab'), printm 'Ab'
	f.nblock = 8;
	Ab = Gblock(A, f.nblock);
end

% check 3c vs analytical
if 0, printm '3l proj check'
	cpu etic
	pp = Ab * xtrue;
	cpu etoc '3l proj'
%	im clf, im_toggle(proj(:,:,1:12:end), pp(:,:,1:12:end), [0 4.4])
	nrms(pp, proj)
end

if 0 % tests
	tmp = Ab{1} * xtrue(ig.mask);
	tmp = Ab{1} * xtrue;
	ia = 1:f.nblock:cg.na;
	tmp = Ab{1}' * col(li_hat(:,:,ia));
	tmp = A * ig.unitv;
	im(tmp), cbar
	tmp = A' * tmp;
	im(tmp), cbar
return
end

% regularization object
if ~isvar('R'), printm 'regularizer'
	f.l2b = 2^4.5;
	f.delta = 100/1000;
	R = Reg1(ig.mask, 'type_denom', 'matlab', ...
                'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b);
	if 1 % check spatial resolution (away from edges)
		W = diag_sp(yi(:));
		psf = qpwls_psf(A, R, 1, ig.mask, W, 'fwhmtype', 'profile');
	end
end

if ~isvar('xpwls'), printm 'PWLS reconstruction'
	si = log(bi ./ yi); % log sinogram
	wi = yi;
	xpwls = pwls_pcg1(xfdk(ig.mask), A, diag_sp(wi), si(:), R, ...
		'niter', 100, 'stop_threshold', 1e-3);
	xpwls = ig.embed(xpwls);
	im(xpwls), cbar
prompt
end

% reshape data to be "2d arrays" for OS iterations (subset over last dim)
if ~isvar('os_data'), printm 'os_data'
	if isscalar(bi) && isscalar(ri)
		os_data = reshaper(yi, '2d');
		os_data = {os_data, ...
			bi * ones(size(os_data)), ri * ones(size(os_data))};
	else
		os_data = {reshaper(yi, '2d'), reshaper(bi, '2d'), ...
			reshaper(ri, '2d')}; % all data as 2d arrays
	end
end

% OS-SPS iterations for transmission penalized likelihood
if ~isvar('xpl'), printm 'start iterations'
	f.niter = 20;
	xinit = xfdk;
	xs = tpl_os_sps(xinit(ig.mask), Ab, os_data{:}, R, 1+f.niter);
	xs = ig.embed(xs);
	xpl = xs(:,:,:,end);
	im(xpl)
end

% finally, compare FDK vs iterative results
if 1
	clim = [0 0.02];
	im clf
	im pl 2 2
	im(1, xtrue, 'true', clim), cbar
	im(2, xfdk, 'FDK', clim), cbar
	im(3, xpl, 'PL', clim), cbar
	im(4, xpwls, 'PWLS', clim), cbar
	nrms(xtrue, xfdk)
	nrms(xtrue, xpl)
	nrms(xtrue, xpwls)
end
