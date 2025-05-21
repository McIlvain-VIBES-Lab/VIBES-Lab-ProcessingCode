function T1_dcm2nii
%% Anatomical Scan DCM to NII
% March 26th 2025
% Grace McIlvain

% This code converts a T1 anatomical to a nii

   %!/Applications/MRIcron/dcm2niix * -4fn
   % !/Applications/MRIcron/dcm2niix *
   % !gunzip -f *.nii.gz
   % !mv co*.nii coT1W_3D_TFE.nii
   % !cp coT1W_3D_TFE.nii ../coT1W_3D_TFE.nii

   %TS added May 2 Note: Some Subjects have 2 T1W and Horos does not handle
   %this well 
    items = dir('T1W_3D_TFE');
    for k = 1:length(items)
        itemName = items(k).name;
        fullPath = fullfile('T1W_3D_TFE', itemName);
        
        if items(k).isdir && ~strcmp(itemName, '.') && ~strcmp(itemName, '..')
            rmdir(fullPath, 's'); 
        end
    end
   !/Applications/MRIcron/dcm2niix -z n -b n T1W_3D_TFE
   !mv T1W_3D_TFE/*.nii coT1W_3D_TFE.nii

end

  
