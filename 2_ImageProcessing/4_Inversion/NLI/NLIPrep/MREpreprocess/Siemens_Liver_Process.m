clear all

% magdir='EPIMRE_928_3D_PRISMA_32SLC_RR_MAG_0029';
% phsdir='EPIMRE_928_3D_PRISMA_32SLC_RR_P_P_0030';
% 
% nph=4;
% nsl=32;
% ndir=3;

magdir='EPIMRE_928_3D_PRISMA_RR_MAG_0017';
phsdir='EPIMRE_928_3D_PRISMA_RR_P_P_0018';

nph=4;
nsl=8;
ndir=3;

  


if(exist('RawData.mat','file'))
    load RawData.mat
else
    % Process the magnitude images
    cd(magdir)
    d=dir('*.IMA');
    for ii=1:length(d)
        info = dicominfo(d(ii).name);
        Y = dicomread(info);
        %figure, imshow(Y);    
        disp(['mag ' int2str(ii) ' ' num2str(info.SliceLocation)])    
        stack(:,:,ii)=Y;    
        SliceLoc(ii)=info.SliceLocation;
    end

    voxsize_mm(1:2)=info.PixelSpacing;
    voxsize_mm(3)=info.SpacingBetweenSlices;

    slicenum=round((SliceLoc-min(SliceLoc))/voxsize_mm(3)+1);
    stack5d=zeros([size(Y) nsl nph ndir]);
    for ii=1:(length(d)/3)
        stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,1)=stack(:,:,ii);
        stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,2)=stack(:,:,ii+nph*nsl);    
        stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,3)=stack(:,:,ii+2*nph*nsl);
    end
    AnatomicalMRI=mean(mean(stack5d,5),4);
    montagestack(AnatomicalMRI)
    title('AnatomicalMRI')

    % Process Phase Images
    cd ..
    cd(phsdir)

    d=dir('*.IMA');
    for ii=1:length(d)
        info = dicominfo(d(ii).name);
        Y = dicomread(info);
        %figure, imshow(Y);    
        disp(['phs ' int2str(ii) ' ' num2str(info.SliceLocation)])    
        stack(:,:,ii)=Y;    
        SliceLoc(ii)=info.SliceLocation;
    end

    slicenum=round((SliceLoc-min(SliceLoc))/voxsize_mm(3)+1);
    stack5d=zeros([size(Y) nsl nph ndir]);
    for ii=1:(length(d)/3)
        stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,1)=stack(:,:,ii);
        stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,2)=stack(:,:,ii+nph*nsl);    
        stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,3)=stack(:,:,ii+2*nph*nsl);
    end

    phs=stack5d/4095*2*pi;

    for ii=1:nph
        for jj=1:ndir
            montagestack(phs(:,:,:,ii,jj))
            title(['Phase ' int2str(ii) ' dir ' int2str(jj)])
            colorbar
        end
    end

    cd ..
    save RawData AnatomicalMRI phs nsl ndir nph voxsize_mm
end

% Masking
if(~exist('Mask.mat','file'))
    mask=Stack_Segmentation(AnatomicalMRI);
    save Mask mask
else
    load Mask.mat
end
image=permute(phs,[1 2 3 5 4]);

% Unwrapping
if(exist('mreimages_unwrap.mat','file'))
    load mreimages_unwrap.mat
else
    [ny,nx,nz,ni,nj] = size(image);
    nk = ni*nj;
    image(isnan(image))=0;
    imrs = reshape(image,[ny nx nz nk]);

    mask_rep = repmat(mask,[1 1 1 nk]);
    mask_nii = make_nii(mask_rep);
    save_nii(mask_nii,'mre_mask.nii')

    phs_nii = make_nii(mask_rep.*imrs);
    save_nii(phs_nii,'mre_phs.nii')

    mag_rep = repmat(AnatomicalMRI,[1 1 1 nk]);
    mag_nii = make_nii(mag_rep);
    save_nii(mag_nii,'mre_mag.nii')

    disp('running prelude in fsl... ')
    tic
    !echo $FSLDIR/bin/prelude
    !ls $FSLDIR/bin/prelude

    !$FSLDIR/bin/prelude -a mre_mag.nii -p mre_phs.nii -o mre_output.nii -m mre_mask.nii -v
    !gunzip -f mre_output.nii.gz
    toc

    tmp = load_nii('mre_output.nii');
    pimg = double(reshape(tmp.img,[ny nx nz ni nj]));
    save mreimages_unwrap.mat pimg ny nx nz ni nj
end
%disp_img = flip(pimg,5)*1.576; % scaling for um, unimportant
fft_img = fft(pimg,[],5)/(nj/2)*1.576; % scaling for um, unimportant;
wave_img = fft_img(:,:,:,:,2);
t2stack=AnatomicalMRI;

Xmotion = wave_img(:,:,:,1);
Ymotion = wave_img(:,:,:,2);
Zmotion = wave_img(:,:,:,3);

mreParams.subj = 'mre_for_inversion';
mreParams.FOVx = nx*voxsize_mm(1);
mreParams.FOVy = ny*voxsize_mm(2);
mreParams.FOVz = nz*voxsize_mm(3);
mreParams.nx = nx;
mreParams.ny = ny;
mreParams.nz = nz;
mreParams.freq = 60;

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack')



