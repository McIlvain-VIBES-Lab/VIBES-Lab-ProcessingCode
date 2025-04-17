% Load the DTI data from WashU and create the DTI.mat file that
% MRE_preprocess needs to run Aniso inversion

% Im assuming that the 3 dimensions of the DTI is in +LPH convention. So
% use the smae conversions as the WashU_to_NLI_convert code

% The DTI and MRE data are different sizes:

% define the position of the DTI (1,1,1) corner in the MRE stack
DTI_111=[3 5 1];

load TBPhantom_20201211_DTI_data_for_Matt.mat
load TBPhantom_20201211_100Hz_LRAct_spiralosc_data_for_Matt.mat

sizDTI=size(evec);
sizMRE=size(PCWU_BMR);


FA1=FA;

FA=zeros(sizMRE(1:3));
V1=zeros([sizMRE(1:3) 3]);

% Need to zero pad if the MRE space doesnt fit inside the DTI space. Will Do this later if required. 
lim1=[DTI_111(1) DTI_111(1)+sizMRE(1)-1];
lim2=[DTI_111(2) DTI_111(2)+sizMRE(2)-1];
lim3=[DTI_111(3) DTI_111(3)+sizMRE(3)-1];
FA=FA1(lim1(1):lim1(2),lim2(1):lim2(2),lim3(1):lim3(2));
V1=evec(lim1(1):lim1(2),lim2(1):lim2(2),lim3(1):lim3(2),[2 1 3]); % Note the direction switch here

%FA1(max(lim1(1),1):min(lim1(2),sizDTI(1)),max(lim1(1),1):min(lim1(2),sizDTI(1)),max(lim1(1),1):min(lim1(2),sizDTI(1)))
montagestack(FA1);
title(['DTI FA size ' int2str(sizDTI(1:3))])

montagestack(MAGmean);
title(['MRE MAGmean size ' int2str(sizMRE(1:3))])

montagestack(FA);
title(['FA and MRE space size ' int2str(sizMRE(1:3))])


montagestack(V1(:,:,:,1));
title(['V1 component 1'])
colorbar

montagestack(V1(:,:,:,2));
title(['V1 component 2'])
colorbar

montagestack(V1(:,:,:,3));
title(['V1 component 3'])
colorbar


%%
quivsl=20;
figure;
imagesc(MAGmean(:,:,quivsl));colormap gray
hold on
quiver3(MAGmean(:,:,quivsl),V1(:,:,quivsl,2),V1(:,:,quivsl,1),V1(:,:,quivsl,3),3)
%quiver3(MAGmean(:,:,quivsl),abs(V1(:,:,quivsl,2)),0*V1(:,:,quivsl,1),0*V1(:,:,quivsl,3),3);  % This gives +y in my MRE coordnicate system
%quiver3(MAGmean(:,:,quivsl),0*(V1(:,:,quivsl,2)),abs(V1(:,:,quivsl,1)),0*V1(:,:,quivsl,3),3) % This gives +x in my MRE coordnicate system. so i think ive got it right.

axis equal

save DTI V1 FA
