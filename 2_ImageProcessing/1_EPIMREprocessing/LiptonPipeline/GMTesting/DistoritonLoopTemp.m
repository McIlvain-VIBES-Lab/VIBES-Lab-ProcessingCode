dirlist = dir('2023*')
for jj=15%56:length(dirlist)
    cd(dirlist(jj).name)
    clear t2stack
    clear magimg
    clear mag
    clear t2bet
    clear Xmotion
    clear Ymotion
    clear Zmotion
    clear mask
    
    t2tmp1 = load_nii('t2stack.nii');
    dx = 1.875;
    dy = 1.875;
    dz = 2.5;
    freq = 50;
    [~, SubjectName] = system('basename "$PWD"');

    SubjectName = strtrim(SubjectName);
    addpath(SubjectName)
clear magimg 
clear phsimg
    cd('File_Storage/Ax_BRAIN_MRE')
    dirlist2 = dir('IM*');
    nn = 0;
    for bb =1:4
        for aa = 1:(length(dirlist2)/4) %% number of slices
            %    for dd = 1:45
            nn=nn+1;
            images1(:,:,aa,bb) = dicomread(dirlist2(nn).name);
        end
    end

    tt = 1;
    rr = 1;
    for jj=1:(length(dirlist2)/48)
        for ii=1:4

            e1_r(:,:,rr,ii) = images1(:,:,tt,ii);
            e1_i(:,:,rr,ii) = images1(:,:,tt+1,ii);

            e2_r(:,:,rr,ii) = images1(:,:,tt+2,ii);
            e2_i(:,:,rr,ii) = images1(:,:,tt+3,ii);

            e3_r(:,:,rr,ii) = images1(:,:,tt+4,ii);
            e3_i(:,:,rr,ii) = images1(:,:,tt+5,ii);

            e4_r(:,:,rr,ii) = images1(:,:,tt+6,ii);
            e4_i(:,:,rr,ii) = images1(:,:,tt+7,ii);

            e5_r(:,:,rr,ii) = images1(:,:,tt+8,ii);
            e5_i(:,:,rr,ii) = images1(:,:,tt+9,ii);

            e6_r(:,:,rr,ii) = images1(:,:,tt+10,ii);
            e6_i(:,:,rr,ii) = images1(:,:,tt+11,ii);

        end
        tt = tt+12;
        rr =rr+1;
    end


    e1_t = double(e1_r)+i*double(e1_i);
    e2_t = double(e2_r)+i*double(e2_i);
    e3_t = double(e3_r)+i*double(e3_i);
    e4_t = double(e4_r)+i*double(e4_i);
    e5_t = double(e5_r)+i*double(e5_i);
    e6_t = double(e6_r)+i*double(e6_i);

    e1 = flip(flip(permute(e1_t,[2 1 3 4]),2),1);
    e2 = flip(flip(permute(e2_t,[2 1 3 4]),2),1);
    e3 = flip(flip(permute(e3_t,[2 1 3 4]),2),1);
    e4 = flip(flip(permute(e4_t,[2 1 3 4]),2),1);
    e5 = flip(flip(permute(e5_t,[2 1 3 4]),2),1);
    e6 = flip(flip(permute(e6_t,[2 1 3 4]),2),1);

    cd ..
    cd('DistortionCorrectionFiles')

    phsimg(:,:,:,3,:)=angle(e1./e2);
    magimg(:,:,:,3,:)=((abs(e1)+abs(e2))/max(col(abs(e1)+abs(e2))))*30000;
    phsimg(:,:,:,2,:)=angle(e3./e4);
    magimg(:,:,:,2,:)=((abs(e3)+abs(e4))/max(col(abs(e3)+abs(e4))))*30000;
    phsimg(:,:,:,1,:)=angle(e5./e6);
    magimg(:,:,:,1,:)=((abs(e5)+abs(e6))/max(col(abs(e5)+abs(e6))))*30000;

    % GM build in distortion correction
    %GE_MRE_DistortionCorrection

    mag = reshape(magimg,[128 128 size(magimg,3) 12]);
    phs = reshape(phsimg,[128 128 size(magimg,3) 12]);

    cplx_img = mag.*exp(1i*phs);

    realimg = real(cplx_img);
    imagimg = imag(cplx_img);

    realimg2 = (flip(permute(realimg,[2 1 3 4]),2));
    imagimg2 = (flip(permute(imagimg,[2 1 3 4]),2));

    real_nii = make_nii(realimg2,[dx dy dz]);
    save_nii(real_nii,'real.nii')

    imag_nii = make_nii(imagimg2,[dx dy dz]);
    save_nii(imag_nii,'imag.nii')

    imgraw = reshape(cplx_img,[128 128 size(cplx_img,3) 1 3 4]);
    save imgraw_ep2d.mat imgraw

    !$FSLDIR/bin/fugue -i real.nii --dwell=.000356 --loadfmap=QSM/FM2MRE.nii --unwarpdir=x- -u QSM/fugue_output_real
    !$FSLDIR/bin/fugue -i imag.nii --dwell=.000356 --loadfmap=QSM/FM2MRE.nii --unwarpdir=x- -u QSM/fugue_output_imag

    cd QSM/
    !gunzip -f *.nii.gz
    cd ../

    real_nii_fix = load_untouch_nii('QSM/fugue_output_real.nii');
    imag_nii_fix = load_untouch_nii('QSM/fugue_output_imag.nii');
    realimg_fix = real_nii_fix.img;
    imagimg_fix = imag_nii_fix.img;

    real1_fix = permute((flip(realimg_fix,2)),[2 1 3 4]);
    imag1_fix = permute((flip(imagimg_fix,2)),[2 1 3 4]);


    real1_fix = double(real1_fix);
    imag1_fix = double(imag1_fix);
    cplx_img2 = real1_fix + i.*imag1_fix;

    realimg = real(cplx_img2);
    imagimg = imag(cplx_img2);
    imgraw2 = reshape(cplx_img2,[128 128 size(cplx_img2,3) 1 3 4]);
    save imgraw2 imgraw2
    %end
magimg2 = abs(squeeze(imgraw2));
magimg2(isnan(magimg2))=0;
  t2stack = mean(mean(magimg2,5),4);
t2nii = make_nii(t2stack,[dy dx dz]);
save_nii(t2nii,'t2stack.nii')
!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.25 -w 1.1

    !gunzip -f t2bet.nii_mask.nii.gz
    !gunzip -f t2bet.nii.gz
    !cp t2bet.nii_mask.nii t2mask.nii
    !rm t2bet_mask.nii

    tmp = load_nii('t2mask.nii');
    seD1 = strel('diamond',0);
    %imerode gets ride of the extra one bad voxel at the edge
    mask = imerode(double(permute(flip(flip(tmp.img,1),2),[2 1 3])),seD1);
    save t2mask_bet.mat mask

    t2stack = flip(flip(permute(t2stack,[2 1 3]),1),2);
    save t2stack.mat t2stack

    phsimg=angle(squeeze(imgraw2));


    % This code makes phase into wave maps for NLI
    load ../maskx.mat

    mask = mask.*abs(1-maskx);
    save t2mask_final.mat mask

    image = flip(flip(permute(phsimg,[2 1 3 4 5]),2),1);
    image(:,:,:,2,:) = (-1)*image(:,:,:,2,:);

    [ny,nx,nz,ni,nj] = size(image);
    nk = ni*nj;
    image(isnan(image))=0;
    imrs = reshape(image,[ny nx nz nk]);

    mask_rep = repmat(mask,[1 1 1 nk]);
    mask_nii = make_nii(mask_rep);
    save_nii(mask_nii,'mre_mask.nii')

    phs_nii = make_nii(mask_rep.*imrs);
    save_nii(phs_nii,'mre_phs.nii')

    mag_rep = repmat(t2stack,[1 1 1 nk]);
    mag_nii = make_nii(mag_rep);
    save_nii(mag_nii,'mre_mag.nii')

    disp('running prelude in fsl... ')
    tic
    !$FSLDIR/bin/prelude -a mre_mag.nii -p mre_phs.nii -o mre_output.nii -m mre_mask.nii -v
    !gunzip -f mre_output.nii.gz
    toc

    tmp = load_nii('mre_output.nii');
    pimg = double(reshape(tmp.img,[ny nx nz ni nj]));
    save mreimages_unwrap.mat pimg
    %OSS_SNR = oss_snr_filter(pimg,[dx dy dz],freq,mask);
    %!cp dcfiles/OSS_SNR.mat ./
    %disp_img = flip(pimg,5)*1.576; %try this to remove hotness maybe
    disp_img = pimg*1.576; % scaling for um, unimportant
    fft_img = fft(disp_img,[],5)/(nj/2);
    wave_img = fft_img(:,:,:,:,2);

    Zmotion = wave_img(:,:,:,1);
    Ymotion = wave_img(:,:,:,2);
    Xmotion = wave_img(:,:,:,3);

    mreParams.subj = 'mre_for_inversion';
    mreParams.FOVx = nx*dx;
    mreParams.FOVy = ny*dy;
    mreParams.FOVz = nz*dz;
    mreParams.nx = ny;
    mreParams.ny = ny;
    mreParams.nz = nz;
    mreParams.freq = freq;

    save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack')

    SubjectName_DC = sprintf('%s_DC',SubjectName);

    load('mre_for_inversion.mat')
    save(sprintf('%s.mat',SubjectName_DC),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack')
    UIUC_data_convert_mcilvain(SubjectName_DC)
    cd(sprintf('%s',SubjectName_DC))
    MRE_preprocess_v9_mcilvain('default',SubjectName_DC)
    eval(sprintf('!mv %s.mat %s/',SubjectName_DC,SubjectName_DC))
    cd ..
    addpath(SubjectName_DC)
    system(sprintf('scp -r %s gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/',SubjectName_DC));
    pause(20)
    insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName_DC, '/hex/', SubjectName_DC, '_voxelmesh'];
    system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))
   % clear all
    cd ..
    cd ..
    cd ..
end
