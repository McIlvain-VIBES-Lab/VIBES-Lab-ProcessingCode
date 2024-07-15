function mre_3Drecon_epilepsy(datadir)

cd(datadir)

if ~exist('SE_FM','dir')
    mkdir SE_FM
    senfile = dir('*fieldmap*');
    if length(senfile) == 1
        eval(sprintf('!mv %s SE_FM/',senfile(1).name))
    else
        disp('Error: no senmap found')
    end
else
    senfile1 = dir('SE_FM/*fieldmap*');
    if length(senfile1) == 0
        senfile = dir('*fieldmap*');
        if length(senfile) == 1
            eval(sprintf('!mv %s SE_FM/',senfile(1).name))
        else
            disp('Error: no senmap found')
        end
    end
end

MRE_recon_3D_2015_07_24_INSIGHT(1);

cd ..