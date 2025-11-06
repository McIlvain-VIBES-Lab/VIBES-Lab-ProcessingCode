function HelenLiptonPrepareForNLI(mreParams, mask, Zmotion, Ymotion, Xmotion, t2stack, OSS_SNR, mask_label)

%% Section 3A — Prepare & Send to NLI (per mask)

% Identify subject name
[~, SubjectNameBase] = system('basename "$PWD"');
SubjectNameBase = strtrim(SubjectNameBase);

% Append mask label so each mask version is unique
SubjectName = sprintf('%s_%s', SubjectNameBase, mask_label);
disp(['Preparing NLI data for Subject: ', SubjectName]);

% Save variables to a new .mat for this mask version
save(sprintf('%s.mat', SubjectName), 'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR');

% Convert to NLI-compatible data format
UIUC_data_convert_mcilvain(SubjectName);

% Run preprocessing
cd(SubjectName);
MRE_preprocess_v9_mcilvain('default', SubjectName);
eval(sprintf('!mv %s.mat %s/', SubjectName, SubjectName));
cd ..

% Transfer to NLI (Grace/Teah/Hailin path setup)
system(sprintf('scp -r %s ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName));
pause(20);

% Define path for submission
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];

% Submit to SLURM via both users (in case)
system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"', insomniapath));
system(sprintf('ssh ts3641@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"', insomniapath));

disp(['✅ Sent mask ', mask_label, ' to NLI successfully.']);
end
