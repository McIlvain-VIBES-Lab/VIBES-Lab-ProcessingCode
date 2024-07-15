% Gcone_test.m
% test the Gcone object
% Copyright 2008-1-1, Jeff Fessler, University of Michigan

pn = jf_protected_names;

%ptypes = {'sf1', 'nn1', 'pd1'};
ptypes = {'sf1', 'sf2', 'sf3'};
if 1 && exist('dd_ge1_mex') == 3
	ptypes{end+1} = 'dd1'; % UM only
	ptypes{end+1} = 'dd2'; % UM only
end
nn = length(ptypes);

% small systems for basic tests
if 1 || ~isvar('A1'), printm 'setup small'
	f.down = 16;
	igs = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 1, 'dz', 0.5, ...
		'offset_x', 2.9, 'offset_y', 3.5, 'offset_z', -3.4, ...
		'mask', 'all-but-edge', ... % trick: needed for hct2
		'down', f.down);

	dfs_list = [0 inf inf]; % arc flat parallel
	dsd_list = [949.075 949.075 inf]; % arc flat parallel
	for kk = 1:length(dfs_list)
		cgs = ct_geom('ge1', 'nt', 320, ...
			'source_z0', 40, 'pitch', 0.5, ... % test helix
			'dfs', dfs_list(kk), ... % arc or flat
			'dsd', dsd_list(kk), ... % fan or parallel beam
			'down', f.down);
		if im, cgs.plot(igs); end

		clear A1 Ac Ah
		for ii=1:nn
			ptype = ptypes{ii};
			if (streq(ptype, 'dd1') || streq(ptype, 'dd2')) ...
				&& isinf(cgs.dsd)
				A1{ii} = [];
				Ac{ii} = [];
				continue
			end
			A1{ii} = Gcone(cgs, igs, 'type', ptype, 'nthread', 1);
			Ac{ii} = Gcone(cgs, igs, 'type', ptype);
			Ah{ii} = Gcone(cgs, igs, 'type', ptype, 'use_hct2', 1);
		end

		for ii=1:nn
			printm('testing type %s dfs=%g dsd=%g', ...
				ptypes{ii}, cgs.dfs, cgs.dsd)
			if isempty(A1{ii}), continue, end
			tester_tomo2(A1{ii}, igs.mask, 'G2', Ac{ii}) % paces
			test_adjoint(A1{ii}, 'big', 1, 'tol', 5e-5)
			test_adjoint(Ac{ii}, 'big', 1, 'tol', 5e-5)
			if pn.has_hct2
				if streq(ptypes{ii}, 'nn1') || streq(ptypes{ii}, 'pd1')
					thresh = 3e-2; % big because of rounding in nn1
				else
					thresh = 7e-6;
				end
				if 1 % mex vs hct2
					xs = igs.mask;
					t1 = A1{ii} * xs;
					t2 = Ah{ii} * xs;
					equivs(t1, t2, 'thresh', thresh)
					b1 = A1{ii}' * t1;
					b2 = Ah{ii}' * t1;
					equivs(b1, b2, 'thresh', thresh)
				end
				if 0 % block
					B1 = Gblock(A1{ii}, 2);
					Bh = Gblock(Ah{ii}, 2);
					t1 = B1{2} * xs;
					t2 = Bh{2} * xs;
					equivs(t1, t2, 'thresh', thresh)
				end
				if 0
					t0 = A1{1} * xs;
					max_percent_diff(t0, t1)
					max_percent_diff(t0, t2)
					pr nrms(t1(:), t0(:))
					pr nrms(t2(:), t0(:))
				end
				tester_tomo2(Ah{ii}, igs.mask, ...
					'equiv_thresh', thresh) % paces
				test_adjoint(Ah{ii}, 'big', 1, 'tol', 5e-5)
			end
		end
	end
end

if ~isvar('x0'), printm 'x0 big'
	f.down = 4;
	igb = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 1, 'dz', 0.5, ...
		'offset_x', 12.9, 'offset_y', 3.5, 'offset_z', -3.4, ...
		'down', f.down);
	ell = [3*igb.dx 5*igb.dx -2*igb.dz ...
		igb.dx*igb.nx/3 igb.dy*igb.ny/4 igb.zfov/4 ...
		0 0 10];
	x0 = ellipsoid_im(igb, ell, 'oversample', 2);
end

% big systems for accuracy tests
if ~isvar('Ab'), printm 'setup big'
	cgb = ct_geom('ge1', 'nt', 320, ...
       		'source_z0', -40, 'pitch', 0.5, ... % test helix
		'down', f.down);
%		'dfs', inf, ... % flat detector
%		'dsd', inf, 'dfs', inf, ... % parallel beam
%		'na', 1, 'orbit_start', 13, ...

	clear Ab Ah
	for ii=1:nn
		ptype = ptypes{ii};
		if streq(ptype, 'dd', 2) && isinf(cgb.dsd)
			Ab{ii} = [];
			continue
		end
		Ab{ii} = Gcone(cgb, igb, 'type', ptype);
	end
end

if ~isvar('ya'), printm 'analytical projections'
	ya = ellipsoid_proj(cgb, ell, 'oversample', 2);
	im clf, im(ya)
prompt
end

if ~isvar('yb'), printm 'discrete projections'
	nrmse = @(x,y) norm(y(:)-x(:)) / norm(x(:)) * 100;
	for ii=1:nn
		if ~isempty(Ab{ii})
			cpu etic
			yb{ii} = Ab{ii} * x0;
			f.time(ii) = cpu('etoc');
			printm('nrmse %s %g %%', ptypes{ii}, nrmse(ya, yb{ii}))
		else
			f.time(ii) = 0;
			yb{ii} = cgb.zeros;
		end
	end
end

if 0, printm 'look at error in worst views'
	im('pl', 2, nn)
	for ii=1:length(ptypes)
		err = yb{ii} - ya;
		tmp = reshape(err, [], cgb.na);
		tmp = sum(tmp.^2); % error in each view
		ia = imax(tmp); % worst view
		im(ii, err(:,:,ia)), cbar h
		titlef('%s ia=%d', ptypes{ii}, ia)
		im('subplot', ii+nn)
		plot(tmp), axis tight
	end
return
end

if 1, printm 'projection profiles'
	it = cgb.nt;
	it = round(cgb.nt/2); % middle
	it = it + [-2 0 2];
	ia = imin(abs(cgb.ad - 45));
%	ia = ceil(ia/2);
	pro = @(y) col(y(:,it,ia));
	arg = [];
	for ii=1:length(ptypes)
		arg = [arg pro(yb{ii})];
	end
	if im
		clf, plot([arg pro(ya)])
		text(10, 200, sprintf('ia=%d', ia))
		text(10, 400, sprintf('ang=%g', cgb.ad(ia)))
		legend(ptypes{:}, 'true')
		axisy(0, 1.2 * max(ya(:)))
		grid
	end
return
end


if 0 % dd1 vs dd2 - they match well
	i_dd1 = strmatch('dd1', ptypes);
	i_dd2 = strmatch('dd2', ptypes);
	im clf, im(yb{i_dd1} - yb{i_dd2}), cbar
	equivs(yb{i_dd1}, yb{i_dd2}, 'thresh', 2e-5)
return
end

if 1, printm 'show projections and differences'
	im clf, im('pl', 2, 1+nn)
	im(1, x0)
	ia = round([1 cgb.na/4 cgb.na/2 cgb.na]);
	im(nn+2, ya(:,:,ia))

	for ii=1:nn
		tmp = yb{ii};
		im(ii+1, tmp(:,:,ia))
		xlabel(ptypes{ii})
		im(ii+2+nn, tmp(:,:,ia)-ya(:,:,ia))
	end
prompt
end

if 1, printm 'show back-projections'
	im clf, im('pl', 1, nn)
	iz = round([1 igb.nz/4 igb.nz/2 igb.nz]);
	iz = 1:2:igb.nz;
	for ii=1:nn
%		tmp = cgb.ones;
		tmp = cgb.zeros; tmp(:,:,20) = 1;
		tmp = Ab{ii}' * tmp;
		im(ii, tmp(:,:,iz))
		xlabel(ptypes{ii})
	end
end

if 0
	im pl 1 3
	ia = round([1 cg.na/4 cg.na/2 cg.na]);
	im row 4
	im(1, ya(:,:,ia));
	im(2, yc(:,:,ia));
	im(3, yc(:,:,ia)-ya(:,:,ia));
	im reset
end
