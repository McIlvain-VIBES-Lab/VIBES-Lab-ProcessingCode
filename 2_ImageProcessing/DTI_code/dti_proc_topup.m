%% 
%% 
function dti_proc_topup

voxel_size = [2 2 2];

dti = load_untouch_nii('Diff_topup.nii');
dtinii = dti.img;
1
save_untouch_nii(dti,'dti_data.nii')
2
% dti.img = mean(dtinii(:,:,:,1:2),4); (Old version)
% % b0nii = make_nii(mean(dtinii(:,:,:,1:2),4),voxel_size);
% save_untouch_nii(dti,'dti_b0.nii')
dti.img = mean(dtinii(:,:,:,1:2),4);
dti.hdr.dime.dim(5) = 1;
% b0nii = make_nii(mean(dtinii(:,:,:,1:2),4),voxel_size);
save_untouch_nii(dti,'dti_b0.nii')

%-f option sets fractional intensity threshold. Default .5, smaller value gives larger brain outline
%-m generate binary brain dti_mask
%-n dont generate the default brain image output

!$FSLDIR/bin/bet2 dti_b0.nii dti -f 0.3 -w 1 -m -n -v
!gunzip -f dti_mask.nii.gz

% PROCESS DTI TENSOR
% note: I left the "bvecs" and "bvals" files in the code folder... so
% change the hard paths on those to match your setup
!$FSLDIR/bin/dtifit -k dti_data.nii -o dti -m dti_mask.nii -r Diff.bvec -b Diff.bval -V
!gunzip -f dti_FA.nii.gz
!gunzip -f dti_MD.nii.gz
!gunzip -f dti_MO.nii.gz
!gunzip -f dti_S0.nii.gz

!gunzip -f dti_L1.nii.gz
!gunzip -f dti_L2.nii.gz
!gunzip -f dti_L3.nii.gz

% creating an AD image (same as L1)
!cp dti_L1.nii dti_AD.nii
3
% creating an RD image (L1+L2)/2
tmp2 = load_untouch_nii('dti_L2.nii');
tmp3 = load_untouch_nii('dti_L3.nii');
tmpRD = tmp2;
tmpRD.img = (tmp2.img+tmp3.img)/2;
4
save_untouch_nii(tmpRD,'dti_RD.nii')

% note: you don't really need the vectors for analysis, but probably for
% displaying the color FA in FSLview
!gunzip -f dti_V1.nii.gz
!gunzip -f dti_V2.nii.gz
!gunzip -f dti_V3.nii.gz
5
dtinew = load_untouch_nii('dti_v1.nii');
dtin = dtinew.img;
6
FA = load_untouch_nii('dti_FA.nii');

FAmask = FA.img>0.25;
dtimg = dtin.*repmat(FAmask,[1 1 1 3]);

FA = FA.img;
save FA.mat FA

dtix = dtimg(:,:,:,1);
dtiy = dtimg(:,:,:,2);
dtiz = dtimg(:,:,:,3);

cd .. 
cd .. 
save dtix.mat dtix
save dtiy.mat dtiy
save dtiz.mat dtiz

dti1 = cat(4,dtix,dtiy);
dtiall = cat(4,dti1,dtiz);
save dtiall.mat dtiall