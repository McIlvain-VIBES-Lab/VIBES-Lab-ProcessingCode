
% MRE_plotv7: Plots Results for MREv7 reconstructions.
% Versions: 
% v1 : Original during initial development of MREv7 - very basic.
% v2 : Interpolates all values back to MR voxels for display.
% v3 : Works with elemental basis support.
% v4 : Major update for reconstuctions using MRE-Zone-v7.04+ reconstruction
%      and MRIhexmeshinterpv1 meshing step. 
%       - Records reconstuction and mesh data in the same structure as the 
%         results.
%       - Uses data saved from the meshing step to interpolate material
%         properties back onto the original MR voxels. 





%close all;clear all;clc

function MRE_v9_process_mtrfile


if(nargin<3) % Default is to not show images 
  showimg=false;
end
%% Choose which Reconstructed Dataset to use

d=dir('*01.nod');
for ii=1:length(d)
    disp(['File ' int2str(ii) ' :: ' d(ii).name])
end
n=input('Number of file to use (default = 1)  >> ');
if(isempty(n))
    n=1;
end   
nodf=d(n).name;
nodstm=nodf(1:end-14);


d=dir([nodstm(1:end-1) '*RE.*prop.01.mtr']);
if length(d)~=0
    mfind=false;
    for ii=1:length(d)
        disp(['File ' int2str(ii) ' :: ' d(ii).name])
    end
else
    d=dir([nodstm(1:end-1) '*RE.*prop.01.mf.mtr']);
    mfind=true;
    for ii=1:length(d)
        disp(['File ' int2str(ii) ' :: ' d(ii).name])
    end
end
clear n
n=input('Number of file to use (default = last)  >> ');
if(isempty(n))
    n=length(d);
end
mtrf=d(n).name;

% Call MRE interpolation function
MRE_plotv9_mask_Linear_Function(nodf,mtrf)

end
    



    



    
    








