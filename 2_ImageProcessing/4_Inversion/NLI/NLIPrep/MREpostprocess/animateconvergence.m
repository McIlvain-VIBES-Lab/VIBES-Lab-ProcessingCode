%% Process all mtr files into mat files

nodf='multi_inclusion_silicone_2mm_70Hz_voxelmesh_G3300.v7.3.inv.iso.incomp.visc_SPoff_SF0p0015_CG2p2.mtrmesh.01.nod';
d=dir([nodf(1:20) '*RE.*prop.01.mtr']);

for ii=1:length(d)
    MRE_plotv7p3_mask_Linear_Function(nodf,d(ii).name)
end

%% Load all mat files
m=dir([nodf(1:20) '*ReconProps.mat']);
load(m(1).name);
RSMstack=zeros([size(RealShear) length(m)]);

for ii=1:length(m)
    load(m(ii).name);
    RSMstack(:,:,:,ii)=RealShear;
end

%% Play as movie

load ../../../Mask.mat
figure
for ii=1:size(RSMstack,4)
    montagestack(RSMstack(:,:,:,ii)./1000,[],'clip','n')
    title(['Iteration ' int2str(ii)])
    colorbar
    pause(0.5)
end

%% Record movie

load ../../../Mask.mat
h=figure;
vidObj = VideoWriter('iterationprogress.avi');
    
    vidObj.FrameRate = 10;
    open(vidObj)
for ii=1:size(RSMstack,4)
    montagestack(RSMstack(:,:,6:35,ii)./1000,[6,5],'clip','n')
    title(['Iteration ' int2str(ii)])
    caxis([0 4.5])
    colorbar
    colormap(gca,'jet')
    set(h,'position',[165   120   984   532]);
    set(h,'color','w');
    pause(0.4)
    drawnow()
    pause(0.4)
    if(ii<10)
    currFrame = getframe(h);
    size(currFrame.cdata)
    writeVideo(vidObj,currFrame);
    end
    
    if(ii<20)
    currFrame = getframe(h);
    size(currFrame.cdata)
    writeVideo(vidObj,currFrame);
    end
    
    currFrame = getframe(h);
    size(currFrame.cdata)
    writeVideo(vidObj,currFrame);
    
    
end
close(vidObj);

