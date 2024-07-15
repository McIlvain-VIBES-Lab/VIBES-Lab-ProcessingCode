  function pn = jf_protected_names
%|function pn = jf_protected_names
%|
%| A serious drawback of the matlab language is that it lacks
%| a protected or local namespace.  Every m-file that is in the path
%| is available to all functions (except those in "private" subdirectories).
%| Users who have their own m-files that happen to have the same names as
%| any of the routines in a toolbox like this one will have problems.
%|
%| To try to overcome this limitation, I created this function in late 2009
%| to serve as a repository of simple functions.
%| To use any of these functions, one types something like
%|	pn = jf_protected_names;
%| and then one can call the functions using
%|	out = pn.fun(arg1, arg2, ...);
%|
%| Copyright 2009-11-21, Jeff Fessler, University of Michigan


pn = strum(struct, { ...
	'has_hct2', @jf_has_hct2, '()';
	'hct_arg', @jf_hct_arg, '(cg, ig)';
	'ind2sub', @jf_ind2sub, '(siz, ind)';
	});

end % jf_protected_names()


%
% jf_ind2sub()
% version with a single matrix output, one dimension per column
%
function subs = jf_ind2sub(st, Nd, ind)
ind = ind(:);
switch length(Nd)
case 2
	[subs(:,1) subs(:,2)] = ind2sub(Nd, ind);
case 3
	[subs(:,1) subs(:,2) subs(:,3)] = ind2sub(Nd, ind);
case 4
	[subs(:,1) subs(:,2) subs(:,3) subs(:,4)] = ind2sub(Nd, ind);
otherwise
        fail 'not done'
end
end % jf_ind2sub()


%
% function jf_hct_arg()
% name/value pairs needed for hct command line
%
function str = HIDEjf_hct_arg(cg, ig)

if isinf(cg.dfs), fail('dfs=inf not done'), end

str = [
sprintf(' nx %d ny %d nz %d', ig.nx, ig.ny, ig.nz) ...
sprintf(' dx %g dy %g dz %g', ig.dx, ig.dy, ig.dz) ...
sprintf(' offset_x %g offset_y %g offset_z %g', ...
	ig.offset_x, ig.offset_y, ig.offset_z) ...
sprintf(' ns %d nt %d na %d', cg.ns, cg.nt, cg.na) ...
sprintf(' ds %g dt %g', cg.ds, cg.dt) ...
sprintf(' offset_s %g offset_t %g', cg.offset_s, cg.offset_t) ...
sprintf(' dfs %g dso %g dsd %g', cg.dfs, cg.dso, cg.dsd) ...
sprintf(' orbit %g orbit_start %g pitch %g', ...
	cg.orbit, cg.orbit_start, cg.pitch) ...
];
end % HIDE
