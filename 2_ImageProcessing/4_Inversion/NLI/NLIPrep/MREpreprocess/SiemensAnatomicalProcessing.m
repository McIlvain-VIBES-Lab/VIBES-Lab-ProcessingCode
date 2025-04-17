%clear all; close all; 
function []=SiemensAnatomicalProcessing(fname)


if(nargin==0) % Priompt for folder name
    files = dir; 
    % Get a logical vector that tells which is a directory.
    subFolders = files([files.isdir]);
    subFolders(1:2)=[];
    % Extract only those that are directories.

    for k = 1 : length(subFolders)
        fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
    end

    fnum=input('Dicom Directory: >>');
    fname=subFolders(fnum).name;
end

    
firstfile=dir([fname '/*.IMA']);
firstfile=[fname '/' firstfile(1).name];
dcminfo=dicominfo(firstfile);


% nslices=regexp(DicomInfo(1).Private_0051_100a,'\d*','Match'); 
% nslices=str2num(nslices{1,end});
% ndir=3;
% venc=regexp(DicomInfo(1).Private_0051_1014,'\d*','Match'); 
% venc=str2num(venc{1,end});
% 
% RawMotion=zeros(DicomInfo(1).Rows,DicomInfo(1).Columns,nslices,ndir,DicomInfo(1).CardiacNumberOfImages); 



[dicom_images]=dir([fname '/*.IMA']); 
for kk=1:length(dicom_images)
    dicom_image(kk)=dicominfo([fname '/' dicom_images(kk).name]); %reads in each dicom header
    stack(:,:,kk)=dicomread([fname '/' dicom_images(kk).name]);     
end

montagestack(stack); colorbar
[fname '/' dicom_images(kk).name]
save([fname '.mat'],'stack','dcminfo'); 