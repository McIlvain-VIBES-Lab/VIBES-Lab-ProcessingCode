% FIrst run mricron and use dcm2nii to convert into gz.nii file

%!gunzip -f *nii.gz
nii_f=dir('*.nii');
tmp = load_nii(nii_f.name);

colorFA1=flip(permute(tmp.img,[2 1 4 3]),1);

colorFA2=permute(tmp.img,[2 1 4 3]);

figure
montage(colorFA1)

figure
montage(colorFA2)

save ColorFA colorFA1
