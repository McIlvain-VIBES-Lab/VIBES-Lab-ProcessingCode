%%% issue: working on making .mat to nii

startup_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode';
run(fullfile(startup_path, 'startup_matlab_general.m'));

dirlist = dir;
for i = 1:length(dirlist)
    if startsWith(dirlist(i).name, '.')
        continue;  
    end
    
    if exist(fullfile(dirlist(i).name,'mre_for_inversion.mat'))

        lipton_data_dir = fullfile(dirlist(i).name, 'LiptonDataTransfer');
        
        mkdir(lipton_data_dir);
    
        nli_outputs = fullfile(dirlist(i).name, 'nli_outputs');
        niis = fullfile(dirlist(i).name, 'niis');
        DampingRatio = fullfile(nli_outputs, 'DR.mat');
        Stiffness = fullfile(nli_outputs, 'Mu.mat');
        %mag = fullfile(niis, 'mag.nii');
        png = fullfile(dirlist(i).name, 'nli_outputs', '*.png');
    
        
        if exist(nli_outputs, 'dir')
                copyfile(png, lipton_data_dir);
                if exist(niis, 'dir')
                    copyfile(DampingRatio, lipton_data_dir);
                    copyfile(Stiffness, lipton_data_dir);
                    %copyfile(mag, lipton_data_dir);
                else
                    disp(['niis directory not found for subject: ', dirlist(i).name]);
                    continue;
                end
        else
            disp(['nli_outputs directory not found for subject: ', dirlist(i).name]);
            continue;
        end
    
        Mu_orig = load(Stiffness);
        Mu = Mu_orig; 
        Mu.Mu = single(Mu.Mu);
        Mu.Mu(isnan(Mu.Mu)) = 0;
        save_nii(Mu,fullfile(lipton_data_dir,'Stiffness.nii'));
    
        DR_orig = load(DampingRatio);     
        DR = DR_orig; 
        DR.DR = single(DR.DR);
        DR.DR(isnan(DR.DR)) = 0;
        save_nii(DR,fullfile(lipton_data_dir,'DampingRatio.nii'));
    
        % mag = load_nii(mag);     
        % Magnitude = mag; 
        % Magnitude.img = single(Magnitude.img);
        % Magnitude.img(isnan(Magnitude.img)) = 0;
        % save_nii(Magnitude,fullfile(lipton_data_dir,'Magnitude.nii'));
        % delete(fullfile(lipton_data_dir, 'mag.nii'));
    
        if exist(fullfile(dirlist(i).name,'File_Storage'),'dir')
            if exist(fullfile(dirlist(i).name,'File_Storage/T1W_3D_TFE'),'dir')
                dicoms = dir(fullfile(subj_id,'*.dcm'));
                correct_subj_id = dicominfo(fullfile(dicoms(1).folder,dicoms(1).name)).PatientID;
                movefile(dirlist(i).name, correct_subj_id);
                rmdir(dirlist(i).name, 's');
            else
                disp(['File_Storage/T1W_3D_TFE directory not found for subject: ', dirlist(i).name,' and being used for renaming']);
                continue;
            end
        else
            disp(['File_Storage directory not found for subject: ', dirlist(i).name]);
            continue;
        end
    else  
    end
end