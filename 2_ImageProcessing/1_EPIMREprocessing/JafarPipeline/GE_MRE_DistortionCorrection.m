function GE_MRE_DistortionCorrection
%% Distortion Correction Lipton Data

% !$FSLDIR/bin/fsl_prepare_fieldmap SIEMENS phs/phs.nii mag/mag_brain.nii fmap_rads 2.46
% !gunzip -f *.nii.gz
%     !$FSLDIR/bin/fugue -i 1_ep2d_mre/mag/mag.nii --dwell=.0004733 --loadfmap=1_fugue/fmap_rads.nii --unwarpdir=y- -u 1_fugue/fugue_output
%     !$FSLDIR/bin/fugue -i 1_ep2d_mre/Ep2d_Code_Output/real.nii --dwell=.0004733 --loadfmap=1_fugue/fmap_rads.nii --unwarpdir=y- -u 1_fugue/fugue_output_real
%     !$FSLDIR/bin/fugue -i 1_ep2d_mre/Ep2d_Code_Output/imag.nii --dwell=.0004733 --loadfmap=1_fugue/fmap_rads.nii --unwarpdir=y- -u 1_fugue/fugue_output_imag
%     !gunzip -f *nii.gz

end