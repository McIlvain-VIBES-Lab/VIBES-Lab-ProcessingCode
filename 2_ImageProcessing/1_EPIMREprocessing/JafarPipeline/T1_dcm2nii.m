function T1_dcm2nii
%% Anatomical Scan DCM to NII
% March 26th 2025
% Grace McIlvain

% This code converts a T1 anatomical to a nii

   %!/Applications/MRIcron/dcm2niix * -4fn
   !/Applications/MRIcron/dcm2niix *
   !gunzip -f *.nii.gz
   if exist('co*.nii')
   !mv co*.nii coT1W_3D_TFE.nii
   !cp coT1W_3D_TFE.nii ../coT1W_3D_TFE.nii
   else
   !mv *.nii coT1W_3D_TFE.nii
   !cp coT1W_3D_TFE.nii ../coT1W_3D_TFE.nii
   end
end

  
