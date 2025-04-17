clear all
unwrap=false;

ndir=3;

% Automatically sort the phase and magnitude directries
d=dir; % 
nmag=0;
nphs=0;
for ii=3:length(d); % 1st 2 are . and ..
    if(d(ii).isdir)
        d2=dir([d(ii).name '/*.IMA']);
        if(~isempty(d2))
            info=dicominfo([d(ii).name '/' d2(1).name]);
            if((info.Private_0051_1016(1)=='M')||(info.Private_0051_1016(4)=='M'))
                nmag=nmag+1;
                magdir(nmag).name=d(ii).name;
            elseif((info.Private_0051_1016(1)=='P')||(info.Private_0051_1016(4)=='P'))
                nphs=nphs+1;
                phsdir(nphs).name=d(ii).name;
            else
                disp(['WARNING: Could not assosicate directory ' d(ii).name ' with mag or phase'])
            end
        else
            disp(['directory: ' d(ii).name ' does not contain .IMA files'])
        end
    end
end
disp([int2str(nmag) ' magnitude directories found']);
for ii=1:nmag
    disp([int2str(ii) ': ' magdir(ii).name])
end
disp([int2str(nphs) ' phase directories found']);
for ii=1:nphs
    disp([int2str(ii) ': ' phsdir(ii).name])
end

if(nmag==0||nphs~=3)
    disp('sort out folders')
    error('Did not find at least 1 mag or exactly 3 phase directories in working directory')
end

% nph=16;
% nsl=12;
% ndir=3;

freqHz=1;

if(exist('RawData.mat','file'))
    load RawData.mat
else
    % Process the magnitude images : Only using the first one for now,
    % would be better to use 3 when available.
    cd(magdir(1).name)
    d=dir('*.IMA');
    info = dicominfo(d(1).name);
    nph=info.CardiacNumberOfImages;
    for ii=1:length(d)
        info = dicominfo(d(ii).name);
        Y = dicomread(info);
        %figure, imshow(Y);    
        disp(['mag ' int2str(ii) ' ' num2str(info.SliceLocation)])    
        stack(:,:,ii)=Y;    
        SliceLoc(ii)=info.SliceLocation;
    end
    nsl=length(unique(SliceLoc));
    
    voxsize_mm(1:2)=info.PixelSpacing;
    %voxsize_mm(3)=info.SpacingBetweenSlices;
    voxsize_mm(3)=info.SliceThickness;

    slicenum=round((SliceLoc-min(SliceLoc))/voxsize_mm(3)+1);
    stack5d=zeros([size(Y) nsl nph ndir]);
%     for ii=1:(length(d)/3)
%         stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,1)=stack(:,:,ii);
%         stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,2)=stack(:,:,ii+nph*nsl);    
%         stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,3)=stack(:,:,ii+2*nph*nsl);
%     end
    for ii=1:length(d)
        stack5d(:,:,slicenum(ii),floor((ii-1)/nsl)+1,1)=stack(:,:,ii);
        %stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,2)=stack(:,:,ii+nph*nsl);    
        %stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,3)=stack(:,:,ii+2*nph*nsl);
    end
    AnatomicalMRI=mean(mean(stack5d,5),4);
    montagestack(AnatomicalMRI)
    title('AnatomicalMRI')

    % Process Phase Images
    stack5d=zeros([size(Y) nsl ndir nph]);
    for pp=1:3
        cd ..
        cd(phsdir(pp).name)

        d=dir('*.IMA');
        
        for ii=1:length(d)
            % Get some useful dicom fields
            info = dicominfo(d(ii).name);
            patpos=info.PatientPosition;
            PEdir=info.InPlanePhaseEncodingDirection;
            posdirec=info.Private_0051_1013;
            numphs=info.CardiacNumberOfImages;
            imgorientation=info.ImageOrientationPatient;
            imgpos=info.ImagePositionPatient;
            scantyp=info.Private_0051_100e;
            imgcomments=info.ImageComments;
            Trig(ii)=info.TriggerTime;
            bf=info.Private_0051_1014;in=findstr(bf,'_'); VENC=str2num(bf(2:in(1)-1)); Direc(pp).str=bf(in(end)+1:end);            
            Y = dicomread(info);
            %imagesc(Y);title(['Image ' int2str(ii) ' Sliceloc ' num2str(info.SliceLocation) 'dir ' int2str(pp)]);caxis([1600 2600]);colorbar;pause(0.01)    
            disp(['phs ' int2str(ii) ' ' num2str(info.SliceLocation) ',Dir ' Direc(pp).str ',VENC ' num2str(VENC) ',scantyp ' scantyp])    
            stack(:,:,ii)=Y;    
            SliceLoc(ii)=info.SliceLocation;
        end

        slicenum=round((SliceLoc-min(SliceLoc))/voxsize_mm(3)+1);
        T=range(Trig);
        Trig=Trig-min(Trig);
        phsnum=1+Trig/T*(nph-1);
        Trigtime=unique(Trig);
    %     for ii=1:(length(d)/3)
    %         stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,1)=stack(:,:,ii);
    %         stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,2)=stack(:,:,ii+nph*nsl);    
    %         stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,3)=stack(:,:,ii+2*nph*nsl);
    %     end
        for ii=1:length(d)
            %floor(((1:25)-1)./24+1)
            stack5d(:,:,slicenum(ii),pp,phsnum(ii))=stack(:,:,ii);
            %stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,2)=stack(:,:,ii+nph*nsl);    
            %stack5d(:,:,slicenum(ii),mod(ii-1,4)+1,3)=stack(:,:,ii+2*nph*nsl);
        end

    end

    phs=stack5d/4095*2*VENC-VENC;   
    cd ..
    save RawData AnatomicalMRI phs nsl ndir nph voxsize_mm Direc VENC scantyp imgorientation Trig PEdir patpos imgcomments Trigtime
end

% Masking
if(~exist('Mask.mat','file'))
    % Use bet2 from FSL
    t2nii = make_nii(flipdim(flipdim(permute(AnatomicalMRI,[2 1 3]),1),2),[voxsize_mm(2) voxsize_mm(1) voxsize_mm(3)]);
    save_nii(t2nii,'t2stack.nii')
     % If mask is too 'round' (cuts off top and bottom0, decrease the -w
     % number
    if(nsl>8) % Use old call
        !$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.25 -w 0.9
        !gunzip -f t2bet.nii_mask.nii.gz
        !gunzip -f t2bet.nii.gz
        !cp t2bet.nii_mask.nii t2mask.nii
        !rm t2bet.nii_mask.nii        
        
    else % use linited z FOV option with bet
        !$FSLDIR/bin/bet t2stack.nii t2bet.nii -Z -m -v -f 0.45 -g 0.9
        !gunzip -f t2bet_mask.nii.gz
        !gunzip -f t2bet.nii.gz
        !cp t2bet_mask.nii t2mask.nii
        !rm t2bet_mask.nii
    
    end

    tmp = load_nii('t2mask.nii');
    seD1 = strel('diamond',2);
    mask = imerode(double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3])),seD1);
    %mask = double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3]));
    
    save t2mask_bet.mat mask
    save AnatomicalMRI.mat AnatomicalMRI    
    %mask=Stack_Segmentation(AnatomicalMRI);
    save Mask mask
    montagestack(AnatomicalMRI.*mask)
    title('bet2 mask created')
    montagestack(AnatomicalMRI)
    title('Full anatomical Image')
    
   
    
else
    load Mask.mat
end



for jj=1:ndir
    phsim=squeeze(phs(:,:,:,jj,:)).*repmat(mask,[1 1 1 nph]);
    montagestack(phsim,[nph,nsl],'y');colorbar
    title(['Phase images dir ' int2str(jj) ' ' Direc(jj).str])
    xlabel('slice')
    ylabel('phase')
    colorbar
end

%image=permute(phs,[1 2 3 5 4]);

[nx,ny,nz,ni,nj] = size(phs);
fft_img = fft(phs,[],5)/(nj/2);
wave_img = fft_img(:,:,:,:,2);

dircount=0;
Ur=zeros([nx ny nz ndir]);
Ui=zeros([nx ny nz ndir]);

disp(['Scan type: ' scantyp])
disp(['image orientation ' num2str(imgorientation')])
disp(['Directions ' Direc(1).str ' ' Direc(2).str ' ' Direc(3).str])
if(strcmp(scantyp(1:3),'Cor')) % Coronal
    % 'through' = AP = (1,1,:) = z From phantom, bottom slice stationalry,
    % top slice moving, positve z = P. 
    % 'rl' = (1,:,1) = y, L to right of image, incresing points to L.  
    % 'fh' = (:,1,1) = x , H to top of image, moving down increasing 1st
    % index, toward foot, so needs a -ve sign to fit LPH convention
    % 
    % Coronal definition: H to top of image and L to right

    for ii=1:3
        if(strcmp(Direc(ii).str,'through'))
            dircount=dircount+1;
            disp('found z direction')
            Ur(:,:,:,3)=real(wave_img(:,:,:,ii));
            Ui(:,:,:,3)=imag(wave_img(:,:,:,ii));
        elseif(strcmp(Direc(ii).str,'rl'))
            dircount=dircount+1;
            disp('found y direction')
            Ur(:,:,:,2)=real(wave_img(:,:,:,ii));
            Ui(:,:,:,2)=imag(wave_img(:,:,:,ii));
        elseif(strcmp(Direc(ii).str,'fh'))
            dircount=dircount+1;
            disp('found x direction')
            Ur(:,:,:,1)=-real(wave_img(:,:,:,ii));
            Ui(:,:,:,1)=-imag(wave_img(:,:,:,ii));
        end
    end
    Motion2Image=eye(3);
    RHcoord=false; % Not a RH coordinaate system 

elseif(strcmp(scantyp(1:3),'Tra')) % Axial I think 
    disp('Axial Scan')
    % through = fh = (1,1,:) = z
    % rl = (1,:,1) = y
    % ap = (:,1,1) = x        
    for ii=1:3
        if(strcmp(Direc(ii).str,'through'))
            dircount=dircount+1;
            disp('found z direction')
            Ur(:,:,:,3)=real(wave_img(:,:,:,ii));
            Ui(:,:,:,3)=imag(wave_img(:,:,:,ii));
        elseif(strcmp(Direc(ii).str,'rl'))
            dircount=dircount+1;
            disp('found y direction')
            Ur(:,:,:,2)=real(wave_img(:,:,:,ii));
            Ui(:,:,:,2)=imag(wave_img(:,:,:,ii));
        elseif(strcmp(Direc(ii).str,'ap'))
            dircount=dircount+1;
            disp('found x direction')
            Ur(:,:,:,1)=real(wave_img(:,:,:,ii));
            Ui(:,:,:,1)=imag(wave_img(:,:,:,ii));
        end
    end
    Motion2Image=eye(3);
    RHcoord=true;


else
    disp('UNKNOWN Scan type') 

end 

%From Siemens:
%Positive velocity encoding  is r>l, a>p, and f>h.
%In our dicom header, you should see LPH written in the text.
if(dircount<3)
        disp(['ERROR!!! ' int2str(dircount) ' of 3 directions found - NLI Input file not created '])
else
    save('NLI_MRE_Input','Ur','Ui','AnatomicalMRI','Motion2Image','freqHz','RHcoord','magdir','phsdir','imgcomments','Trigtime','voxsize_mm')
end

montagestack(Ur(:,:,:,1).*mask,[],'y');title('Real X');colorbar
montagestack(Ui(:,:,:,1).*mask,[],'y');title('Imag X');colorbar
montagestack(Ur(:,:,:,2).*mask,[],'y');title('Real Y');colorbar
montagestack(Ui(:,:,:,2).*mask,[],'y');title('Imag Y');colorbar
montagestack(Ur(:,:,:,3).*mask,[],'y');title('Real Z');colorbar
montagestack(Ui(:,:,:,3).*mask,[],'y');title('Imag Z');colorbar
    

