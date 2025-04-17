% Run Elis NITI parameter sensitivity function on simulated data
% Im running this on the high res 1.5mm data used to create the simulated
% data because thats the only data i have properties, DTI and dispalcements
% at the same resolution

% 
clinDTI=true;

% Im running this from:
% dspf='Shell_50Hzmod6_mtr1_xfacexdir.forward.dsp.mat';
% propf='Shell_50Hzmod6_mtr1_xfacexdir.forward.props.mat';
dspf='1702full_1mmprop_DTI_G3000_300_mreAPbcs_het_0p2xTI_K1e6_50Hz.forward.HexdspResults.mat';
propf='1702full_1mmprop_DTI_G3000_300_mreAPbcs_het_0p2xTI_K1e6_50Hz.forward.props.mat';

load(dspf);  % displcemennts
load(propf);  % Properties BC type doesnt matter for this, just make sure its the right mtrfor the disps.
load ../../DTI.mat % DTI fiber directions

if(clinDTI)
    disp('Assuming clincal DTI convention')
    junk=V1(:,:,:,1);
    V1(:,:,:,1)=V1(:,:,:,2);
    V1(:,:,:,2)=junk;
else
    disp('Assuming sim DTI convention')
end
%load ../../Mask.mat % data mask
mask=propstack_RE(:, :, :, 1)>0;

%%
% Create required arrays
niti_prop(:,:,:,1)=propstack_RE(:,:,:,1)+1i*propstack_IM(:,:,:,1); %mu
niti_prop(:,:,:,2)=propstack_RE(:,:,:,2); %phi
niti_prop(:,:,:,3)=propstack_RE(:,:,:,3); %zeta

% Dsps do not have a zero buffer like the props do. Why would I do that? Fool. 
Ur=zeros([size(propstack_RE,1) size(propstack_RE,2) size(propstack_RE,3) 3]);
Ui=zeros([size(propstack_RE,1) size(propstack_RE,2) size(propstack_RE,3) 3]);

Ur(2:end-1,2:end-1,2:end-1,1)=dsp(1).vals; % Re ux
Ur(2:end-1,2:end-1,2:end-1,2)=dsp(3).vals; % Re uy
Ur(2:end-1,2:end-1,2:end-1,3)=dsp(5).vals; % Re uz

Ui(2:end-1,2:end-1,2:end-1,1)=dsp(2).vals; % Im ux
Ui(2:end-1,2:end-1,2:end-1,2)=dsp(4).vals; % Im uy
Ui(2:end-1,2:end-1,2:end-1,3)=dsp(6).vals; % Im uz

[OSS_SNR,Motion_SNR,OSS_SNR_Dist,Motion_SNR_Dist,oss,ons,varargout]=Strain_SNR_and_NITI_from_phase(squeeze(Ur(:,:,:,:,1)+1i*Ui(:,:,:,:,1)),mask,[1.5 1.5 1.5]/1000,0.5,2,3,V1,niti_prop)
%[OSS_SNR,Motion_SNR,OSS_SNR_Dist,Motion_SNR_Dist,oss,ons,varargout]=Strain_SNR_and_NITI_from_phase_bkup(squeeze(Ur(:,:,:,:,1)+1i*Ui(:,:,:,:,1)),mask,[1.5 1.5 1.5]/1000,0.5,2,3,V1,niti_prop)
%%

montagestack(varargout(:,:,:,1));colorbar % Mu sensitivity 
title(['mu sens : ' dspf]);xlabel(['propf: ' propf]);
montagestack(varargout(:,:,:,2));colorbar % phi sensitivity
title(['phi sens : ' dspf]);xlabel(['propf: ' propf]);
montagestack(varargout(:,:,:,3));colorbar % zeta sensitivity
title(['zeta sens : ' dspf]);xlabel(['propf: ' propf]);