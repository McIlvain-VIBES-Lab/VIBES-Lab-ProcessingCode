%% Converting Dicom Files to Nifti Files

dir1 = dir('ep2d_*'); 

for ii = 1:length(dir1) 
    if dir1(ii).isdir==1 
        cd (dir1(ii).name);
        disp(pwd)
        cd ('mag');
     
        !/Applications/MRIcron/dcm2nii * -4fn
        !gunzip -f *.nii.gz
        !mv *.nii mag.nii
        cd ..
        cd ('phs');
     
        !/Applications/MRIcron/dcm2nii * -4fn
        !gunzip -f *.nii.gz
        !mv *.nii phs.nii
        cd .. 
        cd ..
    end
end
