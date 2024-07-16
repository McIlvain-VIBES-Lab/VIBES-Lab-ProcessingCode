function dti_proc_rotatebvecs_topup

% This is the code I wrote to 1) register DTI to MRE, 2) rotate the bvecs,
% and 3) recalculate the DTI params

% should be run from DTI folder, I think, but will need to be careful about
% how you reference the t2bet.nii from the MRE you are registering

% this registers the S0 to the T2 from MRE
% t2bet.nii called here
!$FSLDIR/bin/flirt -in dti_S0.nii -ref t2bet.nii -out dti2mre.nii -omat dti2mre.mat

% 1) now apply this registration to the original diffusion data (i.e. not
% the outputs)
% t2bet.nii called here
% (note I have never used this function, so check it)
% !$FSLDIR/bin/applyxfm4d Diff_topup.nii t2bet.nii Diff_rotated.nii dti2mre.mat
!$FSLDIR/bin/flirt -in Diff_topup.nii -ref t2bet.nii -out Diff_rotated.nii -applyxfm -init dti2mre.mat

% 2) now we are going to suck in the old bvecs, rotate them, and save them
% back out
[B, diffdelim] = importdata('Diff.bvec');
R = load('dti2mre.mat','-ascii');
Rx = R(1:3,1:3);

ndir = size(B,2);
for ii = 1:ndir
    Bx(:,ii) = Rx*B(:,ii);
end

dlmwrite('Diff_rotated.bvec',Bx,diffdelim)

% 3) now i just apply the same stuff you had done before, but on the
% rotated data
!gunzip -f *.nii.gz
dti = load_untouch_nii('Diff_rotated.nii');
dtinii = dti.img;
dti.img = dtinii(:,:,:,1);
dti.hdr.dime.dim(5) = 1;
save_untouch_nii(dti,'dti_rot_b0.nii')

!$FSLDIR/bin/bet2 dti_rot_b0.nii dti -f 0.3 -w 1 -m -n -v
!mv dti_mask.nii.gz dti_rot_mask.nii.gz
!gunzip -f dti_rot_mask.nii.gz

!$FSLDIR/bin/dtifit -k Diff_rotated.nii -o dti_rot -m dti_rot_mask.nii -r Diff_rotated.bvec -b Diff.bval -V
!gunzip -f dti_rot_FA.nii.gz
!gunzip -f dti_rot_MD.nii.gz
!gunzip -f dti_rot_MO.nii.gz
!gunzip -f dti_rot_S0.nii.gz

!gunzip -f dti_rot_L1.nii.gz
!gunzip -f dti_rot_L2.nii.gz
!gunzip -f dti_rot_L3.nii.gz

!gunzip -f dti_rot_V1.nii.gz
!gunzip -f dti_rot_V2.nii.gz
!gunzip -f dti_rot_V3.nii.gz

dtinew = load_untouch_nii('dti_rot_v1.nii');
dtin = dtinew.img;
FA = load_untouch_nii('dti_rot_FA.nii');

FAmask = FA.img>0.25;
dtimg = dtin.*repmat(FAmask,[1 1 1 3]);

FA = FA.img;
save FA.mat FA


