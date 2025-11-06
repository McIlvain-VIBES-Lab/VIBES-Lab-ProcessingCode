% Modified by Helen Liang
% Oct 24th， 2025

function CheckForNLIOutputs_MultiMask
%% Pulls NLI outputs for all mask versions of each subject
% October 2025 - Modified by Helen Liang
% Works with multi-mask naming convention: SUBJECT_MASKLABEL (e.g., G012_AH)

% ==== Path setup ====
code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
common_code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/';
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
addpath(code_path);
addpath(common_code_path);
startup_matlab_general;

% ==== Identify all subject versions ====
dirlist = dir('*_*');  % expects folders like G012_AH, G012_GM, etc.
dirlist = dirlist([dirlist.isdir]);
dirlist = dirlist(~ismember({dirlist.name}, {'.','..'}));

% Group folders by base subject name (everything before first underscore)
subjectBases = unique(cellfun(@(x) strtok(x, '_'), {dirlist.name}, 'UniformOutput', false));

for s = 1:length(subjectBases)
    subjBase = subjectBases{s};
    fprintf('\n===== Checking NLI outputs for Subject %s =====\n', subjBase);

    % Find all mask versions for this subject
    subjVersions = {dirlist(startsWith({dirlist.name}, subjBase)).name};

    % Create parent folder for all outputs (if not exists)
    if ~exist(subjBase, 'dir')
        mkdir(subjBase);
    end
    cd(subjBase);

    % Prepare NLI_Outputs folder
    if ~exist('NLI_Outputs', 'dir')
        mkdir('NLI_Outputs');
    end

    % === Loop through all mask-labeled versions ===
    for v = 1:length(subjVersions)
        subjVer = subjVersions{v};
        fprintf('--- Pulling back NLI outputs for %s ---\n', subjVer);

        % Create local subfolder for this version
        localOutDir = fullfile('NLI_Outputs', subjVer);
        if ~exist(localOutDir, 'dir')
            mkdir(localOutDir);
        end

        % Remote paths on insomnia
        remoteHexPath = sprintf('/insomnia001/depts/mcilvain/users/mcilvain/%s/hex', subjVer);
        invCheckPath = sprintf('%s/%s_voxelmesh/inv', remoteHexPath, subjVer);

        % Check if inversion results exist on NLI
        checkCmd = sprintf('ssh gm3128@insomnia.rcs.columbia.edu "ls %s/*0100.prop* > /dev/null 2>&1"', invCheckPath);
        [status, ~] = system(checkCmd);

        if status ~= 0
            fprintf('⚠️  No completed inversion found for %s (yet).\n', subjVer);
            continue;
        end

        % Pull back results
        fprintf('✅ Found inversion output. Pulling back %s...\n', subjVer);
        system(sprintf('scp -r gm3128@insomnia.rcs.columbia.edu:%s %s', remoteHexPath, localOutDir));

        % === Post-process the inversion results ===
        cd(localOutDir);
        try
            % Run MRE postprocessing on pulled data
            MRE_v9_process_folder;

            % Navigate into voxelmesh/inv
            cd('hex');
            vmDir = dir('*voxelmesh');
            if isempty(vmDir), error('No voxelmesh folder found.'); end
            cd(vmDir(1).name);
            cd('inv');

            invFiles = dir('*avlast08..0100.Re*');
            if isempty(invFiles)
                error('No inversion result files found for %s', subjVer);
            end

            load(invFiles(1).name);

            % Clean NaN/invalid voxels
            Rx = double(abs(RealShear-3300)<0.0001);
            Rxm = abs(1-Rx);
            RealShear = RealShear.*Rxm;
            ImagShear = ImagShear.*Rxm;
            DR = DR.*Rxm;

            ComplexShear = RealShear + (1i * ImagShear);
            AbsShear = sqrt(RealShear.^2 + ImagShear.^2);
            Mu = 2 * (AbsShear.^2) ./ (RealShear + AbsShear);

            save('RealShear.mat','RealShear');
            save('ImagShear.mat','ImagShear');
            save('DR.mat','DR');
            save('ComplexShear.mat','ComplexShear');
            save('AbsShear.mat','AbsShear');
            save('Mu.mat','Mu');

            % Visualization (optional)
            figure; im(Mu(:,:,:)); caxis([0 6000]); colorbar; colormap(gca,stiff_color);
            print('-dpng','-r300',sprintf('Mu_%s', subjVer));

            cd('../../../../'); % back to NLI_Outputs/<subjVer>
            fprintf('✅ Finished processing %s\n', subjVer);

        catch ME
            fprintf('⚠️  Error processing %s: %s\n', subjVer, ME.message);
            cd('../../../../');
        end

        cd('../..'); % back to subject folder
    end

    cd('..'); % back to main directory
end

disp('✅ All available NLI outputs checked and pulled successfully.');
end
