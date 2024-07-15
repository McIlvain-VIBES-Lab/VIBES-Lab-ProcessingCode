function mre_3Drecon_insight(datadir)

curdir = pwd;
cd(datadir)

if (~exist('recon_complete.txt','file'))&&(~exist('recon_skip.txt','file'))

if ~exist('SE_FM','dir')
    mkdir SE_FM
    senfile = dir('*senmap*');
    if length(senfile) == 1
        eval(sprintf('!mv %s SE_FM/',senfile(1).name))
    else
        disp('Error: no senmap found')
    end
else
    senfile1 = dir('SE_FM/*senmap*');
    if length(senfile1) == 0
        senfile = dir('*senmap*');
        if length(senfile) == 1
            eval(sprintf('!mv %s SE_FM/',senfile(1).name))
        else
            disp('Error: no senmap found')
        end
    end
end

tic;
MRE_recon_3D_2015_07_24_INSIGHT(1,12);
T = toc;

fid = fopen('recon_complete.txt','a+');
fprintf(fid,'Completed: %s \n',date);
fprintf(fid,'Total time: %03.2f mins \n',T/60);
fclose(fid);

end

cd(curdir)