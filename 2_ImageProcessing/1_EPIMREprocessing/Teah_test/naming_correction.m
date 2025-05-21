subjs_dir = '/Volumes/McIlvainDrive2/Lipton_Soccer_Study/sandbox';

% Create a log file to record any problem cases
log_file = fullfile(subjs_dir, 'problem_cases.txt');
fid = fopen(log_file, 'w');

dir_list = dir(subjs_dir);

for i = 1:length(dir_list)
    subj_id = dir_list(i).name;

    if dir_list(i).isdir && ~startsWith(subj_id, '.')
        try
            % create'Archive' and 'Test&Dev' directories
            archive_transfer_path = fullfile(subjs_dir, subj_id, 'Archive');
            test_dev_transfer_path = fullfile(subjs_dir, subj_id, 'Test&Dev');

            if ~exist(archive_transfer_path, 'dir')
                mkdir(archive_transfer_path);
            end

            if ~exist(test_dev_transfer_path, 'dir')
                mkdir(test_dev_transfer_path);
            end

            % name changes and Lipton files created if 'mre_for_inversion.mat' file exists (i.e., post-PreNLI)
            if exist(fullfile(subjs_dir, subj_id, 'mre_for_inversion.mat'), 'file')


                ax_path = regexprep(fullfile(subjs_dir, subj_id, 'File_Storage', 'Ax_Brain_MRE'), '[\\/]+', '/');
                study_path = regexprep(fullfile(subjs_dir, subj_id, 'File_Storage', 'study'), '[\\/]+', '/');
                t1_path = regexprep(fullfile(subjs_dir, subj_id, 'File_Storage', 'T1W_3D_TFE'), '[\\/]+', '/');
                fe_path = regexprep(fullfile(subjs_dir, subj_id, 'File_Storage', 'FE_Phs_and_Mag'), '[\\/]+', '/');
                pe_path = regexprep(fullfile(subjs_dir, subj_id, 'File_Storage', 'PE_Phs_and_Mag'), '[\\/]+', '/');
                ss_path = regexprep(fullfile(subjs_dir, subj_id, 'File_Storage', 'SS_Phs_and_Mag'), '[\\/]+', '/');

                % Check for presence of required subfolders
                has_ax_data = exist(ax_path, 'dir');
                has_fe_data = exist(fe_path, 'dir') && exist(pe_path, 'dir') && exist(ss_path, 'dir');
                has_shared_data = exist(study_path, 'dir') && exist(t1_path, 'dir');

                
                if (has_ax_data && has_shared_data) || (has_fe_data && has_shared_data)

                    % renames using the T1 dicom header
                    t1_dicom = dir(fullfile(t1_path, '*.dcm'));
                    if isempty(t1_dicom)
                        fprintf(fid, '%s: No DICOMs found in T1 folder\n', subj_id);
                        continue;
                    end

                    info = dicominfo(fullfile(t1_path, t1_dicom(1).name));
                    true_subj_id = info.PatientID;

                    if ~strcmp(subj_id, true_subj_id)
                        movefile(fullfile(subjs_dir, subj_id), fullfile(subjs_dir, true_subj_id));
                        name_file = fullfile(subjs_dir, true_subj_id, 'Archieve/orginal_name.txt');
                        fid2 = fopen(name_file, 'w');
                        fprintf(fid2, subj_id);
                        fclose(fid2);
                    end

                    subjs_path = fullfile(subjs_dir, true_subj_id);

                    % Create LiptonTransfer directory for transfering post
                    % NLI 
             
                    lipton_transfer_path = fullfile(subjs_path, 'DataTransfer', 'LiptonTransfer');
                    if ~exist(lipton_transfer_path, 'dir')
                        mkdir(lipton_transfer_path);
                    end

                    nli_outputs_path = fullfile(subjs_path, 'nli_outputs');
                    if exist(nli_outputs_path, 'dir')

                        copyfile(fullfile(subjs_path, 't2stack.nii'), fullfile(lipton_transfer_path, "Magnitude.nii"));
                        copyfile(fullfile(nli_outputs_path, '*'), lipton_transfer_path);

                        t2_path = fullfile(lipton_transfer_path, 'Magnitude.nii');
                        t2 = load_nii(t2_path);
                        system(['bash process_lipton_mre.sh ' t2_path]);
                        
                        %load and save mu, dr and t2 
                        mu = load_nii(t2_path);
                        mu_full = load(fullfile(lipton_transfer_path, 'Mu.mat'));
                        Mu = mu_full.Mu;
                        Mu(isnan(Mu)) = 0;
                        Mu = single(Mu);
                        Mu = flip(flip(permute(Mu, [2 1 3]), 1), 2);
                        mu.img = Mu;
                        mu_path = fullfile(lipton_transfer_path, 'Stiffness.nii');
                        save_nii(mu, mu_path);
                        system(['bash process_lipton_mre.sh ' mu_path]);

                        dr = load_nii(t2_path);
                        dr_full = load(fullfile(lipton_transfer_path, 'DR.mat'));
                        DR = dr_full.DR;
                        DR(isnan(DR)) = 0;
                        DR = single(DR);
                        DR = flip(flip(permute(DR, [2 1 3]), 1), 2);
                        dr.img = DR;
                        dr_path = fullfile(lipton_transfer_path, 'DampingRatio.nii');
                        save_nii(dr, dr_path);
                        system(['bash process_lipton_mre.sh ' dr_path]);

                        % Combine and save all data into one 4D NIfTI file
                        MREall = load_nii(t2_path);
                        MREall.img = cat(4, t2.img, Mu, DR);
                        mre_all_path = fullfile(lipton_transfer_path, 'MREall.nii');
                        save_nii(MREall, mre_all_path);
                        system(['bash process_lipton_mre.sh ' mre_all_path]);

                        % Cleanup 
                        delete(fullfile(lipton_transfer_path, '*.mat'));
                        delete(fullfile(lipton_transfer_path, '*.nii'));
                    else
                        fprintf(fid, '%s: Confirm PostNLI is complete\n', subj_id);
                    end
                else
                    fprintf(fid, '%s: File_Storage missing or has incorrectly named folders for subject\n', subj_id);
                end
            else
                fprintf(fid, '%s: Confirm PreNLI is complete\n', subj_id);
            end
        catch ME
            % log any unexpected error for this subject
            fprintf(fid, '%s: ERROR - %s\n', subj_id, ME.message);
        end
    end
end

% Close the log file
fclose(fid);
