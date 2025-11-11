function JafarCheckForNLIOutputs

code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
common_code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/';
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
addpath(code_path);
addpath(common_code_path)
startup_matlab_general
!rm ._*
dirlist = dir('*Real*');
dirlist = dirlist([dirlist.isdir]);  % keep only directories
dirlist = dirlist(~ismember({dirlist.name}, {'.', '..'}));  % exclude . and ..
for ii=1:length(dirlist)
    
 %   cd(dirlist(ii).name)
    SubjectName = sprintf('%s',dirlist(ii).name);
    
    
 %  if ~exist('NLI_Outputs','dir') %GM 
% remotePath = sprintf('/insomnia001/depts/mcilvain/users/mcilvain/%s/hex', SubjectName);
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];
%system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))

%cd(sprintf('/insomnia001/depts/mcilvain/users/mcilvain/%s/hex/%s_voxelmesh',SubjectName,SubjectName));

filePattern = '*0100.prop*';
checkCmd = sprintf('ssh ts3641@insomnia.rcs.columbia.edu "ls %s/inv/%s > /dev/null 2>&1"', insomniapath, filePattern);
[status, result] = system(checkCmd);

if status == 0
system(sprintf('scp -r ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/%s .', SubjectName));
cd(SubjectName)
!rm ._*
cd hex
!rm ._*
cd(sprintf('%s_voxelmesh',SubjectName));
!rm ._*
cd inv
!rm ._*
cd ../../../
MRE_v9_process_folder

cd hex;
dir2 = dir('*voxelmesh');
cd(dir2.name);
cd inv;
!rm ._*
dir2 = dir('*avlast08..0100.Re*.mat');
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
figure;im(Mu(:,:,:)); caxis([0 4000]); colorbar; colormap(gca,stiff_color);
print('-dpng','-r300',sprintf('Mu_%s',dirlist(ii).name(1:end)))
end

        cd ..


   
%   end
   cd ..
end


% end