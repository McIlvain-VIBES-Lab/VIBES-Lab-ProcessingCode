function dti_proc

voxel_size = [1.6 1.6 1.6];

dti = load_nii('dti.nii');
dtinii = dti.img;

save_nii(dti,'dti_data.nii')

b0nii = make_nii(mean(dtinii(:,:,:,1:2),4),voxel_size);
save_nii(b0nii,'dti_b0.nii')

%-f option sets fractional intensity threshold. Default .5, smaller value gives larger brain outline
%-m generate binary brain dti_mask
%-n dont generate the default brain image output

!$FSLDIR/bin/bet2 dti_b0.nii dti -f 0.3 -w 1 -m -n -v
!gunzip -f dti_mask.nii.gz

% PROCESS DTI TENSOR
% note: I left the "bvecs" and "bvals" files in the code folder... so
% change the hard paths on those to match your setup
!$FSLDIR/bin/dtifit -k dti_data.nii -o dti -m dti_mask.nii -r /Volumes/CLJ-001/Group/drsmitty/ScanData/Multi/Hockey_pilot/Cljohnson_Hockey_Mre/Diff/dti.bvec -b /Volumes/CLJ-001/Group/drsmitty/ScanData/Multi/Hockey_pilot/Cljohnson_Hockey_Mre/Diff/dti.bval -V
!gunzip -f dti_FA.nii.gz
!gunzip -f dti_MD.nii.gz
!gunzip -f dti_MO.nii.gz
!gunzip -f dti_S0.nii.gz

!gunzip -f dti_L1.nii.gz
!gunzip -f dti_L2.nii.gz
!gunzip -f dti_L3.nii.gz

% creating an AD image (same as L1)
!cp dti_L1.nii dti_AD.nii

% creating an RD image (L1+L2)/2
tmp2 = load_nii('dti_L2.nii');
tmp3 = load_nii('dti_L3.nii');
tmpRD = tmp2;
tmpRD.img = (tmp2.img+tmp3.img)/2;
save_nii(tmpRD,'dti_RD.nii')

% note: you don't really need the vectors for analysis, but probably for
% displaying the color FA in FSLview
!gunzip -f dti_V1.nii.gz
!gunzip -f dti_V2.nii.gz
!gunzip -f dti_V3.nii.gz

dtinew = load_nii('dti_v1.nii');
dtin = dtinew.img;

FA = load_nii('dti_FA.nii');
FAmask = FA.img>0.4;
dtimg = dtin.*repmat(FAmask,[1 1 1 3])

dtix = dtimg(:,:,:,1);
dtiy = dtimg(:,:,:,2);
dtiz = dtimg(:,:,:,3);

cd .. 
save dtix.mat dtix
save dtiy.mat dtiy
save dtiz.mat dtiz

dti1 = cat(4,dtix,dtiy);
dtiall = cat(4,dti1,dtiz);
save dtiall.mat dtiall


