
clear all; close all;
delete(gcp('nocreate'));

code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
common_code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/';
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
addpath(code_path);
addpath(common_code_path)
startup_matlab_general

system(['bash /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/Teah_test/nli_push.sh'])

dirlist = dir('*-*');
dirlist = dirlist([dirlist.isdir]);  % keep only directories
dirlist = dirlist(~ismember({dirlist.name}, {'.', '..'})); 

for ii=1:length(dirlist)
    try
        cd(fullfile(dirlist(ii).folder, dirlist(ii).name))
        SubjectName = sprintf('%s',dirlist(ii).name);
    
        needNLI_file = fullfile(dirlist(ii).folder, dirlist(ii).name, 'log', 'needs_postNLI.txt');
        postNLI_file = fullfile(dirlist(ii).folder, dirlist(ii).name, 'log', 'PostNLI_complete.txt');
        if exist(needNLI_file, 'file') && ~exist(postNLI_file, 'file')
            close all;
    
            cd('NLI_Outputs')
                MRE_v9_process_folder
    
                cd hex;
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

                % TS added Oct 9 2025 there is an issue with the versions
                % of Java and Matlab 
                figure('Visible', 'on');           
                set(gcf, 'Renderer', 'painters');  
                
                im(Mu);                          
                caxis([0 6000]);
                colorbar;
                colormap(gca, stiff_color);         
                

                [~, baseName, ~] = fileparts(dirlist(ii).name);
                filename = sprintf('Mu_%s.png', baseName);
                
                % Save the figure safely
                print(gcf, '-dpng', '-r300', filename);

        
                %TS added 56-93 on May 1
                %movefile(SubjectName, 'NLI_Outputs');
                cd ..
                mkdir('Nifti_Data');
        
                load(fullfile('NLI_Outputs', 'Mu.mat'));
                load(fullfile('NLI_Outputs', 'DR.mat'));
                load('t2stack.mat')
        
                Mu(isnan(Mu)) = 0;
                Mu = single(Mu);
                Mu = flip(flip(permute(Mu, [2 1 3]), 1), 2);
        
                DR(isnan(DR)) = 0;
                DR = single(DR);
                DR = flip(flip(permute(DR, [2 1 3]), 1), 2);
        
                t2stack(isnan(t2stack)) = 0;
                t2stack = single(t2stack);
                t2stack = flip(flip(permute(t2stack, [2 1 3]), 1), 2);
        
                Mu_nii = load_nii('t2stack.nii');
                Mu_nii.img = Mu;
                mu_path = fullfile('Nifti_Data', 'Stiffness.nii');
                save_nii(Mu_nii, mu_path);
                system(['bash /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/Teah_test/process_lipton_mre.sh ', mu_path])
        
                DR_nii = load_nii('t2stack.nii');
                DR_nii.img = DR;
                dr_path = fullfile('Nifti_Data', 'DampingRatio.nii');
                save_nii(DR_nii, dr_path);
                system(['bash /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/Teah_test/process_lipton_mre.sh ', dr_path]);
        
                t2_nii = load_nii('t2stack.nii');
                t2_nii.img = t2stack;
                mag_path = fullfile('Nifti_Data', 'Magnitude.nii');
                save_nii(t2_nii, mag_path);
                system(['bash /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/Teah_test/process_lipton_mre.sh ', mag_path]);
        
                delete(fullfile('Nifti_Data', '*.nii'));
    
                logFilePath = fullfile(dirlist(ii).folder,SubjectName,'log/PostNLI_complete.txt');
                fid = fopen(logFilePath, 'a');
                if fid ~= -1
                    fprintf(fid, 'Post NLI Processing Complete for: %s\n', SubjectName);
                    fclose(fid);
                end
                delete(needNLI_file)
        end
        
    catch ME 
        fprintf('Error processing subject %s: %s\n', dirlist(ii).name, ME.message);
        fid_err = fopen(fullfile(dirlist(ii).folder, dirlist(ii).name, 'NLI_Error.txt'), 'a');
        if fid_err ~= -1
            fprintf(fid_err, '[%s] Error processing subject %s:\n%s\n\n', ...
                datestr(now), dirlist(ii).name, getReport(ME, 'extended'));
            fclose(fid_err);
        end
    end
    cd ..
end

        