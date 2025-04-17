%% Compile Lipton Data to Send
% January 5 2025

dirlist = dir('2023*')
for ii=1:length(dirlist)
    cd(dirlist(ii).name)

    load('t2stack.mat')
    load('nli_outputs/DR.mat')
    load('nli_outputs/Mu.mat')
    tmp1 = load_nii('t2stack.nii');
    tmp2 = load_nii('t2stack.nii');
    tmp3 = load_nii('t2stack.nii');

    mkdir('niis')
    tmp1.img = flip(flip(permute(t2stack,[2 1 3]),1),2);
    save_nii(tmp1,'niis/mag.nii')
    tmp2.img = flip(flip(permute(Mu,[2 1 3]),1),2);
    save_nii(tmp2,'niis/Stiffness.nii')
    tmp3.img = flip(flip(permute(DR,[2 1 3]),1),2);
    save_nii(tmp3,'niis/DampingRatio.nii')
    cd ..

    cd('MREforLiptonLab')
    mkdir(dirlist(ii).name)
    cd(dirlist(ii).name)
    tmp1.img = flip(flip(permute(t2stack,[2 1 3]),1),2);
    save_nii(tmp1,'Magnitude.nii')
    tmp2.img = flip(flip(permute(Mu,[2 1 3]),1),2);
    save_nii(tmp2,'Stiffness.nii')
    figure;im(Mu); caxis([0 6000]); colorbar; colormap(gca,stiff_color);  
    print('-dpng','-r300',sprintf('Mu_%s',dirlist(ii).name))
    tmp3.img = flip(flip(permute(DR,[2 1 3]),1),2);
    save_nii(tmp3,'DampingRatio.nii')
   
    cd ../../
end