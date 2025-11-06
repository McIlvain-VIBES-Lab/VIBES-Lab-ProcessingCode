function HelenLiptonCleanData()
%% Section 3B — Clean Up Data (run once after all masks sent)

disp('🧹 Performing final Lipton cleanup after sending all masks...');

mkdir('File_Storage');
!mv Ax_Brain_MRE/ File_Storage
!mv study/ File_Storage
!mv T1W_3D_TFE/ File_Storage/
!mv dcfiles/ File_Storage/
!mv QSM/ File_Storage/

!mv maskx.mat File_Storage/
!mv mre_mag.nii File_Storage/
!mv mre_phs.nii File_Storage/
!mv mre_mask.nii File_Storage/
!mv mre_output.nii File_Storage/
!mv OSS_SNR.mat File_Storage/
!mv mreimages_unwrap.mat File_Storage/
!mv t2mask_bet.mat File_Storage/

disp('✅ Cleanup complete. All raw files archived in File_Storage.');
end
