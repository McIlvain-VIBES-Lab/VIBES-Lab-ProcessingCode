dirs={'topup_mag_phs_00pos';
'topup_real_imag_00pos';
'topup_mag_phs_20pos';
'topup_real_imag_20pos';
'topup_mag_phs_20neg';
'topup_real_imag_20neg';
};
%
%
Ndir = size(dirs,1);
%
nli_suf = '_voxelmesh/inv/';
%nli_suf = '_voxelmesh/inv_SP_1e_09/';
%
Nnli = 100;
Navg = 8;
N_avg_end = 100;
N_avg_start = N_avg_end-Navg+1;
%
startDir = pwd;
for i=1:Ndir
	cd(strcat(dirs{i},'/hex/',dirs{i},nli_suf))
        MREv7_avgRange(N_avg_start,N_avg_end)
        MREoutputMAT((Nnli+1))
        close all
        cd(startDir)
end
