%% NLI Post Processing

%   Use this script to clean up data after you pull it back down from NLI
%   You must be in the folder Date-Caviness to run this

dir1=dir('*.mat'); %you will have to change the name here to the first part of the folder name you created
addpath(pwd)

for ii =1:length(dir1)
    %cd(dir1(ii).name(1:end-4));
    cd hex;
    dir2=dir('*voxelmesh');
    cd(dir2.name);
    cd inv
    %dir2 = dir('*avlast08..0100.Re*');
    dir2 = dir('*');

    for jj =1
        load(dir2(jj).name);
        cd ..
        cd ..
        cd ..
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

        mkdir('nli_outputs')
        save nli_outputs/Mu.mat
        save nli_outputs/DR.mat
        save nli_outputs/Mask.mat    
       
        figure;im(Mu(:,:,:)); caxis([0 6000]); colorbar; colormap(gca,stiff_color);
        print('-dpng','-r300',sprintf('Mu_%s',dir1(ii).name(1:end-4)))
        print('-dpng','-r300',sprintf('nli_outputs/Mu_%s',dir1(ii).name(1:end-4)))


    cd ..
    end
end

% Copy Lipton Data back to directory

dirlist = dir('2023-U7487-0578-MO*.mat')
for ii=1:length(dirlist)
    cd(dirlist(ii).name(1:end-4))
    [~, SubjectName] = system('basename "$PWD"');
    SubjectName = strtrim(SubjectName);
    eval(sprintf('!cp -r nli_outputs /Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA/%s',SubjectName))
    cd ..
end

dirlist = dir('2023-U7487-0578-MO*.mat')
for ii=1:length(dirlist)
    cd(dirlist(ii).name(1:end-4))
    [~, SubjectName] = system('basename "$PWD"');
    SubjectName = strtrim(SubjectName);
    eval(sprintf('!cp -r nli_outputs /Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA/%s',SubjectName))
    cd ..
end