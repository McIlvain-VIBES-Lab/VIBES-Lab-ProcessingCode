
dirlist = dir('2023*')
for ii=26:length(dirlist)
    cd(dirlist(ii).name)
    SubjectName = sprintf('%s_DC',dirlist(ii).name);
    cd('File_Storage/DistortionCorrectionFiles')
   % if ~exist('nli_outputs','dir')
      %  remotePath = sprintf('/insomnia001/depts/mcilvain/users/mcilvain/%s/hex', SubjectName);
        insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];
        %system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))

        filePattern = '*0100.prop*';
        %checkCmd = sprintf('ssh gm3128@insomnia.rcs.columbia.edu "ls %s/inv/%s > /dev/null 2>&1"', insomniapath, filePattern);
        checkCmd = sprintf('ssh gm3128@insomnia.rcs.columbia.edu "ls %s/inv/%s > /dev/null 2>&1"', insomniapath, filePattern);
        [status, result] = system(checkCmd);

    %    if status == 0
            system(sprintf('scp -r gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/%s .', SubjectName));
            cd(SubjectName)
            !rm .DS*
            cd hex
            !rm .DS*
            cd ..
            MRE_v9_process_folder
           
            cd hex;
            !rm .DS*
            dir2 = dir('*voxelmesh');
            cd(dir2.name);
            cd inv;
            dir2 = dir('*avlast08..0100.Re*');
            load(dir2(1).name);
            cd ../../../
        Rx = double(abs(RealShear-3300)<0.0001);
        Rxm = abs(1-Rx);
        RealShear = RealShear.*Rxm;
        ImagShear = ImagShear.*Rxm;
        DR = DR.*Rxm;
        save RealShear.mat RealShear
        save ImagShear.mat ImagShear
        save DR.mat DR

        ComplexShear = RealShear + (i*ImagShear);
        AbsShear = sqrt(((RealShear.^2)+(ImagShear.^2)));
        Mu = 2*(AbsShear.^2)./(RealShear+AbsShear);
        save ComplexShear.mat ComplexShear
        save AbsShear.mat AbsShear
        save Mu.mat Mu
        figure;im(Mu(:,:,:)); caxis([0 6000]); colorbar; colormap(gca,stiff_color);
        print('-dpng','-r300',sprintf('Mu_%s',dirlist(ii).name(1:end)))
        cd ..
        mkdir('nli_outputs')
         save nli_outputs/Mu.mat Mu
         save nli_outputs/DR.mat DR
      %   save nli_outputs/mask.mat mask
         print('-dpng','-r300',sprintf('nli_outputs/Mu_%s',dirlist(ii).name(1:end)))
cd ../../../
        end