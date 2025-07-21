function CheckForNLIOutputs
    code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
    common_code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/';
    addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
    addpath(code_path);
    addpath(common_code_path);
    startup_matlab_general

    dirlist = dir('*_*');
    dirlist = dirlist([dirlist.isdir]);  % keep only directories
    dirlist = dirlist(~ismember({dirlist.name}, {'.', '..'}));  % exclude . and ..

    for ii = 1:length(dirlist)
        cd(dirlist(ii).name)
        SubjectName = dirlist(ii).name;

        if ~exist('NLI_Outputs', 'dir')  % If output folder doesn't exist locally

            % Paths for remote commands
            remoteBasePath = sprintf('/insomnia001/depts/mcilvain/users/mcilvain/%s/hex', SubjectName);
            insomniapath = fullfile(remoteBasePath, [SubjectName, '_voxelmesh']);

            % Submit job remotely
            submitCmd = sprintf(['ssh gm3128@insomnia.rcs.columbia.edu "', ...
                'if [ -d ''%s'' ]; then cd ''%s'' && sbatch McIlvain-Submitv9_visc_incomp; else echo ''Directory does not exist: %s''; fi"'], ...
                insomniapath, insomniapath, insomniapath);
            fprintf('Submitting job with command:\n%s\n', submitCmd);
            [submitStatus, submitOutput] = system(submitCmd);
            fprintf('Job submission status: %d\nJob submission output:\n%s\n', submitStatus, submitOutput);

            % Check if expected output exists remotely
            filePattern = '*0100.prop*';
            checkCmd = sprintf('ssh gm3128@insomnia.rcs.columbia.edu "ls %s/inv/%s > /dev/null 2>&1"', insomniapath, filePattern);
            [status, ~] = system(checkCmd);

            if status == 0
                % Copy remote data locally
                scpCmd = sprintf('scp -r gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/%s .', SubjectName);
                fprintf('Running: %s\n', scpCmd);
                system(scpCmd);
                % The folder 'Ax_BRAIN_MRE' is now in pwd, so cd into it once
                cd(SubjectName);
                MRE_v9_process_folder;


                cd hex;
                dir2 = dir('*voxelmesh');
                if isempty(dir2)
                    error('No voxelmesh directory found inside hex');
                end
                cd(dir2(1).name);
                cd inv;
                dir2 = dir('*avlast08..0100.Re*');
                if isempty(dir2)
                    error('No *avlast08..0100.Re* file found in inv');
                end
                load(dir2(1).name);
                cd ../../../

                % Remove zero/invalid regions
                Rx = double(abs(RealShear - 3300) < 1e-4);
                Rxm = abs(1 - Rx);
                RealShear = RealShear .* Rxm;
                ImagShear = ImagShear .* Rxm;
                DR = DR .* Rxm;

                save('RealShear.mat', 'RealShear');
                save('ImagShear.mat', 'ImagShear');
                save('DR.mat', 'DR');

                ComplexShear = RealShear + 1i * ImagShear;
                AbsShear = sqrt(RealShear.^2 + ImagShear.^2);
                Mu = 2 * (AbsShear.^2) ./ (RealShear + AbsShear);

                save('ComplexShear.mat', 'ComplexShear');
                save('AbsShear.mat', 'AbsShear');
                save('Mu.mat', 'Mu');

                figure; im(Mu(:,:,:));
                caxis([0 7000]);
                colorbar;
                colormap(gca, stiff_color);
                print('-dpng', '-r300', sprintf('Mu_%s', SubjectName));
            else
                fprintf('%s NLI data not found.\n', SubjectName);
            end
        else
            fprintf('%s NLI_Ouputs directory already exists locally. Skipping.\n', SubjectName);
        end

        cd ..  % Return to parent dir
    end
end

