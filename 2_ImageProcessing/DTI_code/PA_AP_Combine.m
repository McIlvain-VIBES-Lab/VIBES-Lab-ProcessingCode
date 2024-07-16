type=input('mac or pc?\n','s');
if strcmp(type,'mac')
    subjects=dir('/Volumes/fsergi-001/OsiriX/MORTON_STROKE_MOTOR_LEARNING/Converted_Data/S*');
    path = '/Volumes/fsergi-001/OsiriX/MORTON_STROKE_MOTOR_LEARNING/Converted_Data/';
else
    subjects=dir('Z:\OsiriX\MORTON_STROKE_MOTOR_LEARNING\Converted_Data\S*');
    path = 'Z:\OsiriX\MORTON_STROKE_MOTOR_LEARNING\Converted_Data\'; 
end

homepath= pwd;

for i = 1:size(subjects,1)
    cd([path subjects(i).name])
    !$FSLDIR/bin/fslmerge -t AP_merged ./*_cmrr_mbep2d_diff_AP_modified_2016*/*.hdr 
    !$FSLDIR/bin/fslmerge -t PA_merged ./*_cmrr_mbep2d_diff_PA_modified_2016*/*.hdr 
    cd(homepath)
end