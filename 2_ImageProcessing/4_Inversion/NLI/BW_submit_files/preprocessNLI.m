startDir = pwd;
%
dirs={'topup_mag_phs_00pos';
'topup_real_imag_00pos';
'topup_mag_phs_20pos';
'topup_real_imag_20pos';
'topup_mag_phs_20neg';
'topup_real_imag_20neg';
};
%
Ndir = size(dirs,1);
%
for i=1:Ndir
	cd(dirs{i})
	MRIhexmesh_v7p3_matt('default')
	cd(startDir)
end
cd ..
