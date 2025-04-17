function [OSS_SNR,Motion_SNR,OSS_SNR_Dist,Motion_SNR_Dist,oss,ons]=Strain_SNR(filtwidth,ord,Nfit)

%% Octahedral Shear Strain Based SNR calcualtion
% Calculates Ocathedral Shear and normal strains (OSS and ONS), the curl of
% the displacements, and an OSS based SNR, along with the motion SNR. These
% SNR measures are based on an average over the whole volume, the
% Distribution of the SNR distribution is also output. Motion SNR's are
% given as an average over the real and imaginary componenets of all three
% directions.
%
% Function Input Arguments (all optional, uses default values if not suppiled):
% filtwidth = width of gaussian filter used for preprocessing the data. default=0.5
% ord = order of fitted polynomial to estiamte strains (default 2, quadratic)
% Nfit = size of local block for polynomial fitting (default 3, i.e. 3x3x3 block of data)% 
%
% Requires Input files Files in the current directory: 
% Either MotionData.mat, HeaderData.mat ErrorMap.mat, Mask.mat 
% OR: Hexahedral mesh files
%
% Addtional Functions Required: montagestack.m, getfilename.m
%
% Will Prompt which type of input to use.
%
% Info when using MRE .mat files
% Noise distribution is generated using ErrorMap.mat.
% As of 29 March 2011 there are 3 different noise estimates. 
% Form 1: The original
% Uses mean(misfit), which I decided was not a good noise measure, because 
% std(misfit) is better. These are stored in a variable ErrorMap 
% in ErrorMap.mat, and does not have the variable stdform. This form is in
% datasets generated before 29 march 2011.
% Form 2: Standard Devation form
% Uses std(misfit) for each voxel. These are stored in a variable ErrorMap,
% in ErrorMap.mat, and are identified by an additional variable stdform=true.
% Note a factor of 0.7979 can be used to convert between form 1 and form 2, 
% try r=randn(10000000,1);mean(abs(r))
% Form 3: Propagation of Errors form. 
% This is my preferred form, and is calcualted using propagation of
% Errors/uncertainties. Briefly, the noise of the original MR images is
% estimated using air regions, this is propagated to noise in the MR phase
% (std(phase)=1/abs(MR) * std(MR)), then this is propagated through the
% Discrete Fourier transform to give noise in the fundamental fourier mode,
% (Motion Amplitude) using std(Re{A})=std(Im{A})=std(abs{A})=sqrt(2/N)*std(phase)
% 
% All these error measures are always positive, so a random sign is applied
% to generate the motion noise distribution, so that the strain noise can 
% be calcualted. It would be better to use another propagation of errors 
% calculation with the OSS formula, but this would be really really 
% difficult, because of the numerical differentiation operators inside the 
% OSS formula. I will attempt it one day. - I did this.
%
% Info when using the Hexahedral mesh files input option:
% It will ask for a disp file, and a noisy disp file. The noise is then 
% calcualted by Noisy Disp - Disp. My idea with this was to use it for 
% simulations, where the original dsp file is from a forward solve, and the
% noisy dsp file is after adding noise. 
%
% Original Author: 
% Matt McGarry
% Thayer School of Engineering, Dartmouth College
% 29 March 2011
%
% Version History:
% Original Version, Strain_SNR.m Matt McGarry, 29 March 2011.
% Derivatives estimated by fitting local blocks of data with a generalized
% 3D polynoimal, and differentating that polynomial. Data is prefiltered
% using a 3d Gaussian, note that this filtering is turned off for the noise
% calcualation. A bunch of options are at the start of the code, check
% these out.
%
% % Update 15 Sept 2011, Matt McGarry: 
% Uses propagation of errors to estimate the OSS noise rather than the
% randomized sign approach used previously. Things to note about this is 
% that it assumes a cubic voxel size, inaccuracies may result with highly 
% anistropic voxels, and the noise estimate assumes a central difference 
% approximation of the strains. This code actually uses a polynomial fit 
% (savitsky-golay filter) for calculating the derivatives required to 
% calculate the octahedral shear strain. Also, noise values are not
% available near the borders of the mask. 
% 


%% Preferences

% Select which dataset to use by 'fold'.
fold=[pwd '/']; % Current folder

Parr= true; %true; % Paralell Matlab indicator
nlabs=2; % Number of Matlab Processes to use

opt=3; % Output Option (1=lots of figs, 2=just OSS, 3=no figures output)

ord_def=2; % Order of polynomial
Nfit_def=3; %Size of fitting block

%% Start of Code

%% Define Input Type:
%intyp=input('Input Data Type <1=MRE .mat files, default>, <2 = Hexahedral Mesh Dsp Files>');
intyp=[];disp('Note: Choice for hex mesh disabled, uncomment line above this one to enable')
if(isempty(intyp))
    intyp=1;
end

% Process Filtering/Fitting Parameters
if(nargin<1) % use defaults for filter
    filtwidth=0.5; % Width of Gaussian smoothing filter (pixels). Wider = more smoothing.     
    disp(['Using Default Filter (width = ' num2str(filtwidth) ')'])
end
if(nargin<2)
    ord=ord_def; % order of polynomial to fit (only 2nd/quadratic 3rd/cubic or 4th/quartic currently supported.
    disp(['Using Default Polynomial order :: ' int2str(ord)]) 
end
if(nargin<3)
    Nfit=Nfit_def; % Size of block to fit data to. Needs to be odd.
    disp(['Using Default Fitting Block Size :: ' int2str(Nfit)])
end

if(mod(Nfit,2)==0)
    error('Nfit must be odd')
end

% Size of gaussian filter kernel (pixels). Needs to be at least 6x as
% wide as filtwidth, and needs to be odd.
filtsiz=ceil(6*filtwidth);
if(filtsiz<3)
    filtsiz=3;
end
if(mod(filtsiz,2)==0) % Filtsiz = even, make it odd.
    filtsiz=filtsiz+1;
end
disp(['Gaussian filter size = ' int2str(filtsiz)])




if(intyp==1) % Use MRE_3DMotionData.mat and ErrorMap.mat
    %% True Displacements
    load([fold 'MRE_3DMotionData.mat']);
    load([fold 'HeaderData.mat']);
    load([fold 'Mask.mat']);
    junk=mask;mask=junk;clear junk % renaming Needed for parallel toolbox to know what type of variable.
    dim=size(mask);
    res=DirIndex(4,1:3).*1e-3;
    
    % Note that I am assuming A is in units of MICRONS. The scaling doenst 
    % really matter, particularly for us as we dont do any absolute value 
    % corrections because it doesnt affect our reconstructions. If 
    % comparisons between absolute strain values from different data sources
    % are to be made, this scaling needs to be addressed. The SNR
    % calcualted is valid no matter what scaling is used.
    A=A.*1e-6;
    Ur=A.*cos(P);
    Ui=A.*sin(P);
        
    % Map MR space measurements to physical xyz coordinate system
    % Uses the transformation dim1 = x, dim2 = y, dim3 = z
    % CONFUSING DIRINDEX PROBLEM
    %  In this code, I assume that:
    %  increasing the 1st index = +ve x-direction
    %  increasing the 2nd index = +ve y-direction
    %  increasing the 3rd index = +ve z-direction, 
    %  which is the same assumption made in the hexahedral meshing code. Eli
    %  and Phil worked out the correct transformation matrix for the hex code
    %  is MPSto123=rowswap*DirIndex(1:3,4:6). If i use the same  transformation 
    %  in this code, the octahedral strains are not right (i.e. do not agree 
    %  with the hex and tet code octahedral strain distributions). The Sinkus 
    %  shear modulus estimate is OK.
    %  If I do not use the rowswap matrix; MPSto123=DirIndex(1:3,4:6), the 
    %  calculated octahedral strains match the hex and tet octahdral strains, 
    %  and the Sinkus shear modulus also looks ok, but different.
    %  So I have gone with no rowswap, because the OSS is correct.

    rowswap=[0 1 0;1 0 0;0 0 1]; % Matrix to swap top 2 rows of a matrix
                                 % This is not required to make this oss = tet
                                 % code oss.
    MPSto123=DirIndex(1:3,4:6); % RSM = OK, oct shear looks good
    %MPSto123=rowswap*DirIndex(1:3,4:6);  % RSM looks good, oct shear does not

    % Perform DirIndex Transformation
    rowUr=zeros(prod(dim),3);
    rowUi=zeros(prod(dim),3);
    for ii=1:3
        rowUr(:,ii)=reshape(Ur(:,:,:,ii),[prod(dim) 1]);
        rowUi(:,ii)=reshape(Ui(:,:,:,ii),[prod(dim) 1]);
    end
    rowUr=rowUr*MPSto123'; % From hexahedral Meshing code
    rowUi=rowUi*MPSto123'; % From hexahedral Meshing code
    % Note the transpose is required because the DirIndex is based on the
    % assumption of premultiplying a 3xNn array of displacements 
    %[ux(1) ... ux(n)]            [u1(1) ... u3(n)]
    %[uy(1) ... uy(n)] = DirIndex*[u2(1) ... u2(n)]   (uxyz=DirIndex*u123)
    %[uz(1) ... uz(n)]            [u3(1) ... u3(n)].
    % If we wish to postmutiply an Nnx3 array, we need to add a tranpose,
    % (uxyz)' = (DirIndex*u123)' = u123'*DirIndex'

    for ii=1:3
        Ur(:,:,:,ii)=reshape(rowUr(:,ii),[dim(1) dim(2) dim(3) 1]);
        Ui(:,:,:,ii)=reshape(rowUi(:,ii),[dim(1) dim(2) dim(3) 1]);
    end
    clear rowUr rowUi
    
    %% Noise Displacements
    
    % Note that 2 different Error estimates are available: 
    % Standard deviation of Misfits (SDM):
    % ErrorMap(i,j,k,dir) = std(Up_dir(i,j,k,:) - real(U(i,j,k,dir)*exp(i*t)))
    % Where Up_dir(i,j,k,:) is the N measurements at each phase offset, and
    % U(i,j,k,dir) is the fitted motion amplitude (assumed to the fundamental
    % mode of the DFT in this case), and t is the location of the phase 
    % offsets, usually 0:2*pi/N:(N-1)*2*pi/N. dir is the direction (1=x,
    % 2=y, 3=z). Note that it needs to be multiplied by sqrt(2/N) to 
    % transform it into noisein the complex motion amplitude.  
    % 
    % Propagation of Errors (PoE)
    % The Noise level is estimated using the MR noise, assuming all motion
    % noise is a result of MR noise. 
    % ErrorMap(i,j,k,dir) = sqrt(2/N) * 1/(|X+iY|) * MRnoise, where X+iY is
    % the magnitude of the complex MR voxel at point ijk.
    %
    % Experiments with phantom data have found SDM and POE are similar in
    % the absensce of motion artifacts, i.e. when motion noise is only 
    % caused by MR noise. When Motion artifacts are present, the error in 
    % the motion amplitude is only captured by the SDM errorMap. Motion 
    % artifacts are important in MRE, any noise measure should be able to 
    % capture their effects, therefore the SDM motion error will be used to 
    % estimate strain SNR. The difference between SDM and PoE is a good 
    % indicator of motion articacts, where SDM >> PoE, artifacts such as 
    % unwrapping errors and Intrra voxel phase dispersion (IVPD) are likely 
    % to be responsible.  
    %
    

    load([fold 'ErrorMap.mat']);

    % Find out what type of errormap this is (mean(misfit), std(misfit), Prop
    % of Errors
    if(exist('stdform','var'))
        if(stdform==true)
            disp('ErrorMap in standard deviation of misfits form')
        else
            disp('ErrorMap in mean of misfits form, applying conversion factor')
            ErrorMap=ErrorMap./0.7979;
        end
    else
        disp('ErrorMap in mean of misfits form, applying conversion factor')
        ErrorMap=ErrorMap./0.7979;
    end
    
    % smooth the SDM error to reduce the variability caused by taking the 
    % standard deviation over a small number of offsets. Experiments showed
    % this moderate smoothing has a relativly small effect on the overall
    % SNR (An increase of 1%). The smoothed distribution will have less
    % effects from the limited sample size for taking the standard deviation.
    
    for ii=1:3
        ErrorMap(:,:,:,ii)=smooth3dimage(ErrorMap(:,:,:,ii),0.5,mask); % use a gaussian smoothing filter with SD=1.5
    end
    
    % Assume N=8, apply factor of sqrt(2/N) for Phase error to Motion
    % amplitude error conversion <- This comes from the propagation of 
    % errors analysis. Errormap is an estimate of std(phasevals). The 
    % amplitude error is sqrt(2/N)*std(phasevals). 
    Ur_n=0.5.*ErrorMap.*1e-6;
    Ui_n=0.5.*ErrorMap.*1e-6;

    % Use the Propagation of Errors form of the ErrorMap if availiable
    if(exist('ErrorMapPropofErrs','var'))
        disp('Propagation of Errors form of Errormap present')
        % Real Imaginary Amplitude Noise is equal, Assume ErrorMap and A
        % are in microns. As long as assumption is consistent the actual
        % scaling factor does not affect the SNR.
        UrPOE_n=ErrorMapPropofErrs.*1e-6;
        UiPOE_n=ErrorMapPropofErrs.*1e-6;
        
        % Build an estimate of the motion artifacts
        POEvsSDM = Ur_n./UrPOE_n.*repmat(mask,[1 1 1 3]);
        thresh_art=2; % value 
        figure(99)
        pos = [6 625 1270 296];
        for ii=1:3
            subplot(1,3,ii)
            montagestack(POEvsSDM(:,:,:,ii),[],'clip','n')
            colorbar
            colormap(jet)
            title(['Direction ' int2str(ii)])
            caxis([1.5 10])
        end
        suptitle('Artifact Level estiamte')
        xlabel('values >~2 indicate a probable motion artifact')
        set(gcf,'position',pos)
        save ArtifactEstimate POEvsSDM 
    else
        disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
        disp('No Propagation of Errors ErrorMap Present')
        disp('artifact estimation not performed')
        disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')        
    end
    
    
    % Perform DirIndex Transformation

    MPSto123=DirIndex(1:3,4:6); 
    
    rowUr=zeros(prod(dim),3);
    rowUi=zeros(prod(dim),3);
    for ii=1:3
        rowUr(:,ii)=reshape(Ur_n(:,:,:,ii),[prod(dim) 1]);
        rowUi(:,ii)=reshape(Ui_n(:,:,:,ii),[prod(dim) 1]);
    end
    rowUr=rowUr*MPSto123'; 
    rowUi=rowUi*MPSto123'; 
    % Note the transpose is required because the DirIndex is based on the
    % assumption of premultiplying a 3xNn array of displacements 
    %[ux(1) ... ux(n)]            [u1(1) ... u3(n)]
    %[uy(1) ... uy(n)] = DirIndex*[u2(1) ... u2(n)]   (uxyz=DirIndex*u123)
    %[uz(1) ... uz(n)]            [u3(1) ... u3(n)].
    % If we wish to postmutiply an Nnx3 array, we need to add a tranpose,
    % (uxyz)' = (DirIndex*u123)' = u123'*DirIndex'

    for ii=1:3
        Ur_n(:,:,:,ii)=reshape(rowUr(:,ii),[dim(1) dim(2) dim(3) 1]);
        Ui_n(:,:,:,ii)=reshape(rowUi(:,ii),[dim(1) dim(2) dim(3) 1]);
    end
    clear rowUr rowUi

else % Use Hexahedral Mesh Files
    
    nodf=getfilename('Node File Name','*.nod');
    idxf=getfilename('Index File Name','*.idx');
    dspf=getfilename('Displacement File Name','*.dsp');
    dspnsef=getfilename('Noisy Displacement File Name','*dsp*');
    
    nod=load(nodf);
    nod=nod(:,2:4);
    idx=load(idxf);
    idx=idx(:,2:4);
    dsp=load(dspf);
    dsp=dsp(:,2:7);
    dspnse=load(dspnsef);
    dspnse=dspnse(:,2:7);
    % Calculate Noise by subtracting noisy displacements from normal 
    dspnse=dspnse-dsp;
    
    % Find Resolution
    res=zeros(1,3);
    for ii=1:3
        v=unique(nod(:,ii));
        v=sort(v);
        res(ii)=mean(diff(v));
    end
    
    
    mxind=max(idx);
    % Generate Ur and Ui, Ur_n and Ui_n, and mask
    Ur=zeros(mxind);
    Ui=zeros(mxind);
    Ur_n=zeros(mxind);
    Ui_n=zeros(mxind);
    mask=zeros(mxind);
    jord=[2 1 3]; % Reverse the effect of 'Rowswap' in the hexahedral meshing code.
                  % This makes the OSS distribution the same as the
                  % tetrahedral OSS calculation, but I am not quite sure
                  % why this is.
    for ii=1:length(idx)
        for jj=1:3
            Ur(idx(ii,1),idx(ii,2),idx(ii,3),jord(jj))=dsp(ii,2*jj-1);
            Ui(idx(ii,1),idx(ii,2),idx(ii,3),jord(jj))=dsp(ii,2*jj);
            Ur_n(idx(ii,1),idx(ii,2),idx(ii,3),jord(jj))=dspnse(ii,2*jj-1);
            Ui_n(idx(ii,1),idx(ii,2),idx(ii,3),jord(jj))=dspnse(ii,2*jj);            
        end
        mask(idx(ii,1),idx(ii,2),idx(ii,3))=1;
    end
end

%% Octahedral Shear Strain Calculations

runtme=tic; %start timer

% if isempty(gcp('nocreate'))
%     parpool
% end

% OSS of motions
[oss,ons]=OSS_calc(Ur,Ui,mask,res,filtsiz,filtwidth,ord,Nfit,Parr);

%% OSS Noise Estiamte 
% Use the function OSS_noise, which generates the OSS noise using a 
% propagation of errors analysis, assuming a central difference 
% approximation of the derivatives. Also assumes isotropic voxels. 

% Check isotropic voxels
h=mean(res);
dh=range(res);
% if(dh/h)>0.05 %% GM commented this out on 3/26/25
%     disp(['XXX Warning: voxel anistropy > 5%, dh/mean(h) = ' dh/h])
%     cont=input('XXX Continue using average voxel dimension for noise estimate <y(default)/n>??','s');
%     if(~isempty(cont))
%         if(strcmp(cont,'n'))
%             error('SNR calcualtion aborted because of anisotropic voxels')
%         end
%     end
% end
% Use Gaussian integration to calcualte the average value of OSS noise.
a=0;
b=pi;
wt(1)=(322-13*sqrt(70))/900;
wt(2)=(322+13*sqrt(70))/900;
wt(3)=128/225;
wt(4)=wt(2); %(322+13*sqrt(70))/900;
wt(5)=wt(1); %(322-13*sqrt(70))/900;

gp(1)=-(1/3*sqrt(5+2*sqrt(10/7)));
gp(2)=-(1/3*sqrt(5-2*sqrt(10/7)));
gp(3)=0.d0;
gp(4)=-gp(2); %(1/3*sqrt(5-2*sqrt(10/7)));
gp(5)=-gp(1); %(1/3*sqrt(5+2*sqrt(10/7)));
xi=(a+b)/2.d0+(b-a)/2.d0*gp;

oss_nse=zeros(size(mask));
ossCD=zeros(size(mask));

cdtic=tic;
if(Parr)
    oss_noisegp=zeros([size(mask) 5]);
    ossCDgp=zeros([size(mask) 5]);
    Ugp=zeros([size(mask) 5]);
    Vgp=zeros([size(mask) 5]);
    Wgp=zeros([size(mask) 5]);    
    noisemaskgp=zeros([size(mask) 5]);
    parfor ii=1:5
        % Calculate motion components at the gauss point
        Ugp(:,:,:,ii)=Ur(:,:,:,1).*cos(xi(ii)) - Ui(:,:,:,1).*sin(xi(ii));
        Vgp(:,:,:,ii)=Ur(:,:,:,2).*cos(xi(ii)) - Ui(:,:,:,2).*sin(xi(ii));
        Wgp(:,:,:,ii)=Ur(:,:,:,3).*cos(xi(ii)) - Ui(:,:,:,3).*sin(xi(ii));

        % Std deviation of motion noise is constant over the harmonic cycle
        [oss_noisegp(:,:,:,ii),ossCDgp(:,:,:,ii),noisemaskgp(:,:,:,ii)]=OSS_noise(Ugp(:,:,:,ii),Vgp(:,:,:,ii),Wgp(:,:,:,ii),h,Ur_n(:,:,:,1),Ur_n(:,:,:,2),Ur_n(:,:,:,3),mask);
    end
    for ii=1:5
        oss_nse=oss_nse+oss_noisegp(:,:,:,ii)*wt(ii)*((b-a)/2)/pi;
        ossCD=ossCD+ossCDgp(:,:,:,ii)*wt(ii)*((b-a)/2)/pi; % Use this to check the central diff scheme against the savitsky-golay scheme.
    end
    noisemask=noisemaskgp(:,:,:,1);
else
    for ii=1:5
        % Calculate motion components at the gauss point
        Ugp=Ur(:,:,:,1).*cos(xi(ii)) - Ui(:,:,:,1).*sin(xi(ii));
        Vgp=Ur(:,:,:,2).*cos(xi(ii)) - Ui(:,:,:,2).*sin(xi(ii));
        Wgp=Ur(:,:,:,3).*cos(xi(ii)) - Ui(:,:,:,3).*sin(xi(ii));

        % Std deviation of motion noise is constant over the harmonic cycle
        [oss_noisegp,ossCDgp,noisemask]=OSS_noise(Ugp,Vgp,Wgp,h,Ur_n(:,:,:,1),Ur_n(:,:,:,2),Ur_n(:,:,:,3),mask);
        oss_nse=oss_nse+oss_noisegp*wt(ii)*((b-a)/2)/pi;
        ossCD=ossCD+ossCDgp*wt(ii)*((b-a)/2)/pi; % Use this to check the central diff scheme against the savitsky-golay scheme.
    end
end
disp(['CD oss time = ' num2str(toc(cdtic))])



%% OSS SNR
Imsk=find(noisemask==1);
OSS_SNR=mean(oss(Imsk))/mean(oss_nse(Imsk));
OSS_SNR_Dist=oss./oss_nse;
OSS_SNR_Dist(noisemask~=1)=0;

Motion_SNR=0;
Motion_SNR_Dist=zeros(size(mask));
for ii=1:3
    junk=Ur(:,:,:,ii);
    Sig=mean(abs(junk(Imsk)));   
    junk=Ur_n(:,:,:,ii);
    Nse=mean(abs(junk(Imsk)));    
    Motion_SNR=Motion_SNR+1/6*Sig/Nse;
    Motion_SNR_Dist=Motion_SNR_Dist + 1/6.*abs(Ur(:,:,:,ii))./abs(Ur_n(:,:,:,ii));
    
    junk=Ui(:,:,:,ii);
    Sig=mean(abs(junk(Imsk)));   
    junk=Ui_n(:,:,:,ii);
    Nse=mean(abs(junk(Imsk)));    
    Motion_SNR=Motion_SNR+1/6*Sig/Nse;
    Motion_SNR_Dist=Motion_SNR_Dist + 1/6.*abs(Ui(:,:,:,ii))./abs(Ui_n(:,:,:,ii));
end
OSS_SNR_Dist(mask~=1)=0;
%% Call Output Subroutine
%  Note: output format specifier, opt set at top of code.
disp('Outputting Figures. Change the variable opt at the beginning')
disp('of the code to output figures in a different format')
outputfigs(Ur,Ui,oss,ons,fold,opt)
disp(' ')

    
disp(['Reconstruction walltime = ' num2str(toc(runtme)) ' seconds'])       


% Check out new OSS noise estimate and such like
if(0==1)
    figure(100)
    montagestack(oss,[],'n','n')
    colorbar
    title('oss - SG filtering')
    figure(101)
    montagestack(ossCD,[],'n','n')
    colorbar
    title('oss - central differencing')
    figure(102)
    montagestack(oss_nse,[],'n','n')
    colorbar
    title('oss noise')
    save Junk_OSS oss oss_nse ossCD
    
end
 

save OSS_SNR.mat filtwidth Nfit ord OSS_SNR Motion_SNR OSS_SNR_Dist Motion_SNR_Dist oss ons

end % End of Strain SNR calculation

function [oss,ons]=OSS_calc(Ur,Ui,mask,res,filtsiz,filtwidth,ord,Nfit,Parr)

dim=size(mask);

%% Apply initial filters to data
% apply narrow gaussian filter to remove high frequency noise
H = gaussian3d(filtsiz,filtwidth);
junk1=tic;
if(Parr)  % Parrelelized filtering   
    parfor ii=1:3
        %Ur(:,:,:,ii)=imfilter(Ur(:,:,:,ii),H,'symmetric');
        %Ui(:,:,:,ii)=imfilter(Ui(:,:,:,ii),H,'symmetric');
        % filtering function which does not include values outside the mask
        Ur(:,:,:,ii)=filter3Dmask(Ur(:,:,:,ii),mask,H);
        Ui(:,:,:,ii)=filter3Dmask(Ui(:,:,:,ii),mask,H);    
    end
else % Not Parrelelized filtering
    for ii=1:3
        %Ur(:,:,:,ii)=imfilter(Ur(:,:,:,ii),H,'symmetric');
        %Ui(:,:,:,ii)=imfilter(Ui(:,:,:,ii),H,'symmetric');
        % filtering function which does not include values outside the mask
        Ur(:,:,:,ii)=filter3Dmask(Ur(:,:,:,ii),mask,H);
        Ui(:,:,:,ii)=filter3Dmask(Ui(:,:,:,ii),mask,H);    
    end
end
    
disp(['Pre-filtering of data complete, ' num2str(toc(junk1)) ' seconds'])

% Add a zero buffer on all sides of the data, big enough for the
% overlap of the block used to fit the data. Make room for one expansion of
% the data fitting block
lf=floor(Nfit/2); %Limits either side of centre
bufsiz=lf+1; % Size of zero buffer
jnk=Ur;
Ur=zeros([dim+2*bufsiz,3]);
for ii=1:3
    [Ur(:,:,:,ii)]=zerobuffer3d(jnk(:,:,:,ii),bufsiz,'add');
end
jnk=Ui;
Ui=zeros([dim+2*bufsiz,3]);
for ii=1:3
    [Ui(:,:,:,ii)]=zerobuffer3d(jnk(:,:,:,ii),bufsiz,'add');
end
clear jnk
[mask]=zerobuffer3d(mask,bufsiz,'add');
oss=zeros(dim+2*bufsiz);
ons=zeros(dim+2*bufsiz);
curldsp=zeros([dim+2*bufsiz 3])+1i*zeros([dim+2*bufsiz 3]);

%% Fit a polynomial to blocks of data, to estimate derivatves
if(ord==2) %Quadratic poly has 10 terms
    ntrms=10;
elseif(ord==3);
    ntrms=20; % Cubic poly has 20 terms
elseif(ord==4);
    ntrms=35; % Quartic poly has 35 terms
end

nmsk=sum(mask(:)); %Total number of points to process
disp('Beginning Octahedral Shear Strain Calculation')
tinv=tic;

if(Parr)
    parfor ii=1+bufsiz:dim(1)+bufsiz   
        % The loop was too complicated for matlab's automatic parallel variable
        % classifiation, putting the loop body into a seperate function
        % solved this problem.
        [oss(ii,:,:),ons(ii,:,:),curldsp(ii,:,:,:)]=process_slice(ii,Ur,Ui,mask,ord,lf,dim,bufsiz,ntrms,res);       
    end
else
    for ii=1+bufsiz:dim(1)+bufsiz          
        [oss(ii,:,:),ons(ii,:,:),curldsp(ii,:,:,:)]=process_slice(ii,Ur,Ui,mask,ord,lf,dim,bufsiz,ntrms,res);       
    end

end

disp(['Octahedral Shear Strain calculation complete, ' num2str(toc(tinv)) ' seconds'])
fprintf('\n');% New line

% Remove zero buffers
mask=zerobuffer3d(mask,bufsiz,'remove');
oss=zerobuffer3d(oss,bufsiz,'remove');
ons=zerobuffer3d(ons,bufsiz,'remove');

junk=curldsp;
curldsp=zeros([dim 3])+1i*zeros([dim 3]);
for ii=1:3
    curldsp(:,:,:,ii)=zerobuffer3d(junk(:,:,:,ii),bufsiz,'remove');
end

%% Save Results
%eval(['save ' fold 'OSS_SNR.mat filtwidth filtsiz Nfit ord oss ons curldsp'])

% % Determine if a hex dataset exists
% if(exist([fold 'hex'],'dir'))
%     % Get hex folder names
%     H=dir([fold 'hex']);
%     if(length(H)>3)
%         disp(['Only creating Initial guess file for first mesh : ' H(3).name])
%     end
%     if(H(3).isdir)
%         hfold=H(3).name;
%         F=dir([fold 'hex/' hfold '/*.idx']);
%         if(isempty(F))
%             disp('No .idx files present, No initial guess file created')
%         elseif(length(F)>1)
%             disp([int2str(length(F)) ' idx files present, (need exactly 1), no initial guess file created'])
%         else
%             disp('Creating Initial Guess file')
%             junk=load([fold 'hex/' hfold '/' F(1).name]);
%             idx=junk(:,2:4);
%             Hnn=length(idx);
%             
%             % Viscoelastic initial guess file
%             IGn=zeros(Hnn,4);
%             curln=zeros(Hnn,6);
%             octsn=zeros(Hnn,1);
%             
%             for ii=1:Hnn
%                 IGn(ii,1)=IGr(idx(ii,1),idx(ii,2),idx(ii,3));
%                 IGn(ii,2)=IGi(idx(ii,1),idx(ii,2),idx(ii,3));
%                 IGn(ii,3)=1000;
%                 IGn(ii,4)=-1e-5;
%                 curln(ii,1)=real(curldsp(idx(ii,1),idx(ii,2),idx(ii,3),1));
%                 curln(ii,2)=imag(curldsp(idx(ii,1),idx(ii,2),idx(ii,3),1));
%                 curln(ii,3)=real(curldsp(idx(ii,1),idx(ii,2),idx(ii,3),2));
%                 curln(ii,4)=imag(curldsp(idx(ii,1),idx(ii,2),idx(ii,3),2));
%                 curln(ii,5)=real(curldsp(idx(ii,1),idx(ii,2),idx(ii,3),3));
%                 curln(ii,6)=imag(curldsp(idx(ii,1),idx(ii,2),idx(ii,3),3));
%                 octsn(ii)=oss(idx(ii,1),idx(ii,2),idx(ii,3));
%             end    
%             fid=fopen([fold 'hex/' hfold '/'  F(1).name(1:end-4) '.IGvsc'],'w');
%               fprintf(fid,'%7i %12.4e %12.4e %12.4e %12.4e \n',[(1:Hnn)' IGn]');
%             fclose(fid);
% 
%             % Rayleigh Initial Guess file (Rayleigh comp = 0.5)
%             IGn(:,2)=IGn(:,2)./2;
%             IGn(:,4)=-IGn(:,2)./IGn(:,1).*IGn(:,3);
%             fid=fopen([fold 'hex/' hfold '/' F(1).name(1:end-4) '.IGray'],'w');
%               fprintf(fid,'%7i %12.4e %12.4e %12.4e %12.4e \n',[(1:Hnn)' IGn]');
%             fclose(fid);
%             
%             fid=fopen([fold 'hex/' hfold '/'  F(1).name(1:end-4) '.curl'],'w');
%               fprintf(fid,'%7i %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e \n',[(1:Hnn)' curln]');
%             fclose(fid);
%             
%             fid=fopen([fold 'hex/' hfold '/'  F(1).name(1:end-4) '.oss'],'w');
%               fprintf(fid,'%7i %12.4e\n',[(1:Hnn)' octsn]');
%             fclose(fid);
%         end
%     else
%         disp(['No Initial Guess file created, ' H(3).name ' is not a directory'])
%     end   
% else
%     disp('No hex datasets present')
% end

end


function [oss_ii,ons_ii,curldsp_ii]=process_slice(ii,Ur,Ui,mask,ord,lf,dim,bufsiz,ntrms,res)

%% Process slice ii of the data.
    
dim2=size(Ur);
oss_ii=zeros(dim2(2),dim2(3));
ons_ii=zeros(dim2(2),dim2(3));
curldsp_ii=zeros(dim2(2),dim2(3),3);

ar=zeros(ntrms,3);
ai=zeros(ntrms,3);
msft=zeros(3,2);

for jj=1+bufsiz:dim(2)+bufsiz
    for kk=1+bufsiz:dim(3)+bufsiz
        if mask(ii,jj,kk)==1
            W=mask(ii-lf:ii+lf,jj-lf:jj+lf,kk-lf:kk+lf); % Weight matrix = mask, do not fit to data outside mask.
            
            % Scale derivatives according to Resolution
            dx=res(1); % 
            dy=res(2);
            dz=res(3);

            % Octahedral shear strain calc.
            if(sum(W(:))<=2*ntrms) % Insufficent values inside mask to fit smooth polynomial, increase block size for 
                                   % this point to include more data.
                for ll=1:3
                    W=mask(ii-lf-1:ii+lf+1,jj-lf-1:jj+lf+1,kk-lf-1:kk+lf+1);
                    [ar(:,ll),msft(ll,1)]=weightedfit3dpoly(Ur(ii-lf-1:ii+lf+1,jj-lf-1:jj+lf+1,kk-lf-1:kk+lf+1,ll),ord,W);
                    [ai(:,ll),msft(ll,2)]=weightedfit3dpoly(Ui(ii-lf-1:ii+lf+1,jj-lf-1:jj+lf+1,kk-lf-1:kk+lf+1,ll),ord,W);
                end
            else
                for ll=1:3
                    [ar(:,ll),msft(ll,1)]=weightedfit3dpoly(Ur(ii-lf:ii+lf,jj-lf:jj+lf,kk-lf:kk+lf,ll),ord,W);
                    [ai(:,ll),msft(ll,2)]=weightedfit3dpoly(Ui(ii-lf:ii+lf,jj-lf:jj+lf,kk-lf:kk+lf,ll),ord,W);
                end
            end
            [oss_ii(jj,kk),ons_ii(jj,kk),curldsp_ii(jj,kk,:)]=octshearstain_curl(ar,ai,dx,dy,dz);
        end
    end
end

end

function [ossval,onsval,curlu]=octshearstain_curl(ar,ai,dx,dy,dz)
% Function to calculate the octahedral shear strain given a polynomial fit
% of a local block of dispalcements
% Inputs : ar = Polynomial fit coefficents of local blocks of Real displacements
%          ai = Polynomial fit coefficents of local blocks of Imaginary displacements
% Output : ossval = average octahedral shear strain at center of fitted local block
%          onsval = average octahedral normal strain at center of fitted local block.
%          curlu = curl of displacement field



% % Determine order of polynomial fit from number of coefficents
% Ntrm=length(ar);
% if(Ntrm==4);
%     ord=1;
% elseif(Ntrm==10)
%     ord=2;
% elseif(Ntrm==20)
%     ord=3;
% else
%     error(['Ntrms = ' int2str(length(Ntrm)) ' not recognised'])
% end


% Strain_normal = 1/3 * (S_xx +S_yy + S_zz)
% Strain_shear = 2/3 * [ (S_xx-S_yy)^2 + (S_yy-S_zz)^2 + (S_zz-S_xx)^2 + 3/2*(Sxy^2 + Syz^2 + S_zx^2) ]^(1/2)

% e.g. cubic poly for 5x5x5 data
% f(x,y,z) = a1 + a2*x +a3*y + a4*z + a5*x^2 + a6*y^2 + a7*z^2 
%            + a8*xy + a9*xz + a10*yz + a11*x^3 + a12*y^3 + a13+z^3 
%            + a14*x^2y + a15*xy^2 + a16*x^2z + a17*xz^2 + a18*y^2z + a19*yz^2 + a20*x*y*z
% 1+3+6+10 = 20

dudx=(ar(2,1)+1i*ai(2,1))/dx;
dvdx=(ar(2,2)+1i*ai(2,2))/dx;
dwdx=(ar(2,3)+1i*ai(2,3))/dx;
dudy=(ar(3,1)+1i*ai(3,1))/dy;
dvdy=(ar(3,2)+1i*ai(3,2))/dy;
dwdy=(ar(3,3)+1i*ai(3,3))/dy;
dudz=(ar(4,1)+1i*ai(4,1))/dz;
dvdz=(ar(4,2)+1i*ai(4,2))/dz;
dwdz=(ar(4,3)+1i*ai(4,3))/dz;
sxx=dudx;
syy=dvdy;
szz=dwdz;
sxy=dudy+dvdx;
sxz=dudz+dwdx;
syz=dvdz+dwdy;
% Set up 5 point gaussian integration to get average strain (See Oct Shear
% SNR Paper for why I do this.
a=0;
b=pi;
wt(1)=(322-13*sqrt(70))/900;
wt(2)=(322+13*sqrt(70))/900;
wt(3)=128/225;
wt(4)=wt(2); %(322+13*sqrt(70))/900;
wt(5)=wt(1); %(322-13*sqrt(70))/900;

gp(1)=-(1/3*sqrt(5+2*sqrt(10/7)));
gp(2)=-(1/3*sqrt(5-2*sqrt(10/7)));
gp(3)=0.d0;
gp(4)=-gp(2); %(1/3*sqrt(5-2*sqrt(10/7)));
gp(5)=-gp(1); %(1/3*sqrt(5+2*sqrt(10/7)));

xi=(a+b)/2.d0+(b-a)/2.d0*gp;

octsgp=2/3*sqrt(real((sxx-syy)*exp(1i*xi)).^2+real((syy-szz)*exp(1i*xi)).^2+real((szz-sxx)*exp(1i*xi)).^2 +3/2*(real(sxy*exp(1i*xi)).^2+real(sxz*exp(1i*xi)).^2+real(syz*exp(1i*xi)).^2));
ossval=(sum(wt.*octsgp)*(b-a)/2)/pi; % Average value over one cycle
onsval=1/3*abs(sxx+syy+szz)*2/pi; % Average absolute value over one cycle

curlu(1)=dwdy-dvdz;
curlu(2)=dudz-dwdx;
curlu(3)=dvdx-dudy;

end


function [stackout]=filter3Dmask(stackin,mask,H)
% Performs spatial filtering but does not use values outside the Mask

dim=size(stackin);
stackout=stackin;
filtsiz=size(H);
lf=floor(filtsiz/2); % Distance back and forward of centre of filter
cn=ceil(filtsiz/2); % centre of filter
if(mod(filtsiz(1),2)==0||mod(filtsiz(2),2)==0||mod(filtsiz(3),2)==0)
    error('Filter must be odd-sized')
end

for ii=1:dim(1)
    for jj=1:dim(2)
        for kk=1:dim(3)
            if(mask(ii,jj,kk)==1) %perform filtering
                lm=[max(1,ii-lf(1)) min(dim(1),ii+lf(1)); % Limits to cut out subregion
                      max(1,jj-lf(2)) min(dim(2),jj+lf(2));
                      max(1,kk-lf(3)) min(dim(3),kk+lf(3))];
                subreg=stackin(lm(1,1):lm(1,2),lm(2,1):lm(2,2),lm(3,1):lm(3,2));
                submask=mask(lm(1,1):lm(1,2),lm(2,1):lm(2,2),lm(3,1):lm(3,2));
                % Need to build filter same size as subregion
                la=[max(1,cn(1)-ii+1) max(1,cn(2)-jj+1) max(1,cn(3)-kk+1)]; % first value in each dimension..
                ss=size(subreg);
                Hc=H(la(1):la(1)+ss(1)-1,la(2):la(2)+ss(2)-1,la(3):la(3)+ss(3)-1);
                Hmsk=submask.*Hc;
                Hmsk=Hmsk./sum(Hmsk(:));
                stackout(ii,jj,kk)=sum(Hmsk(:).*subreg(:));                
            end             
        end
    end
end

end


function [arrayout]=zerobuffer3d(arrayin,bufsiz,task)
% Adds or removes a zero buffer around a 3D array
% Inputs: arrayin : 3-D array to be buffered
%         bufsiz : Size of zero buffer to be added to each side of arrayin.
%         task : 'add' will add a buffer, 'remove' will remove it again. (Default = 'add')
% Author: Matt McGarry, 9 June 2010.
%         matthew.d.mcgarry@dartmouth.edu
if(nargin<3)
    task='add';
end

siz=size(arrayin);
if(strcmp(task,'add'))
    newsiz=siz+2*bufsiz;
    arrayout=zeros(newsiz);
    arrayout(bufsiz+1:siz(1)+bufsiz,bufsiz+1:siz(2)+bufsiz,bufsiz+1:siz(3)+bufsiz)=arrayin;
elseif(strcmp(task,'remove'))
    arrayout=arrayin(bufsiz+1:siz(1)-bufsiz,bufsiz+1:siz(2)-bufsiz,bufsiz+1:siz(3)-bufsiz);
end

end

function [a,msft]=weightedfit3dpoly(data,ord,W)
% Fits a 3D polynomial of order=ord to the data submatrix (must be odd-sized),
% Inputs: data(n,n,n) :: 3D block of data to fit 3D polynomial to. (n odd)
%         ord: order of polynomial to fit (4 maximum)
%         W(n,n,n) :: 3D weighting matrix. This can be used to apply no
%                     weight to values outside the mask for polynomial fitting.
% Outputs: a: vector of polynomial coefficents. Length depends on ord.
%          msft: Weight adjusted average misfit of fitted polynomial.


siz=size(data);
if(length(siz))~=3
    error('Requires 3D Data')
end
if~((siz(1)==siz(2))&&(siz(2)==siz(3)))
    error('Requires equal dimension data')
end
if(mod(siz(1),2)==0)
    error('Dimensions of data must be odd')
end
N=siz(1);

% Check Weighting Matrix
if size(W)~=size(data)
    error('Weight matrix must be same dimensions as data')
end
Wm=diag(W(:));

lf=floor(N/2); % number of pixels back and forward of centre

[x,y,z]=meshgrid(-lf:lf,-lf:lf,-lf:lf);
if(ord==1)
    nterm=4;
elseif(ord==2)
    nterm=10;
elseif(ord==3)
    nterm=20;
elseif(ord==4)
    nterm=35;
else
    error('Max supported polynomial order exceeded')
end

a=zeros(nterm,1); %Solution vector of polynomial coefficients
M=zeros(numel(data),nterm);


% e.g. cubic poly for 5x5x5 data
% f(x,y,z) = a1 + a2*x +a3*y + a4*z + a5*x^2 + a6*y^2 + a7*z^2 
%            + a8*xy + a9*xz + a10*yz + a11*x^3 + a12*y^3 + a13+z^3 
%            + a14*x^2y + a15*xy^2 + a16*x^2z + a17*xz^2 + a18*y^2z + a19*yz^2 + a20*x*y*z
% 1+3+6+10 = 20

% Constant term
M(:,1)=1;
% Linear terms
if(ord>=1) 
    M(:,2)=x(:);
    M(:,3)=y(:);
    M(:,4)=z(:);
end
% quadratic terms
if(ord>=2)
    M(:,5)=x(:).*x(:);
    M(:,6)=y(:).*y(:);
    M(:,7)=z(:).*z(:);
    M(:,8)=x(:).*y(:);
    M(:,9)=x(:).*z(:);
    M(:,10)=y(:).*z(:);
end
% cubic terms
if(ord>=3)
    M(:,11)=x(:).*x(:).*x(:);
    M(:,12)=y(:).*y(:).*y(:);
    M(:,13)=z(:).*z(:).*z(:);
    M(:,14)=x(:).*x(:).*y(:);
    M(:,15)=x(:).*y(:).*y(:);
    M(:,16)=x(:).*x(:).*z(:);
    M(:,17)=x(:).*z(:).*z(:);
    M(:,18)=y(:).*y(:).*z(:);
    M(:,19)=y(:).*z(:).*z(:);
    M(:,20)=x(:).*y(:).*z(:);
end
% quartic terms
if(ord>=4)
    M(:,21)=x(:).*x(:).*x(:).*x(:);
    M(:,22)=y(:).*y(:).*y(:).*y(:);
    M(:,23)=z(:).*z(:).*z(:).*z(:);
    M(:,24)=x(:).*x(:).*x(:).*y(:);
    M(:,25)=x(:).*x(:).*x(:).*z(:);
    M(:,26)=y(:).*y(:).*y(:).*x(:);
    M(:,27)=y(:).*y(:).*y(:).*z(:);
    M(:,28)=z(:).*z(:).*z(:).*x(:);
    M(:,29)=z(:).*z(:).*z(:).*y(:);
    
    M(:,30)=x(:).*x(:).*y(:).*y(:);
    M(:,31)=x(:).*x(:).*z(:).*z(:);
    M(:,32)=y(:).*y(:).*z(:).*z(:);
    
    M(:,33)=x(:).*x(:).*y(:).*z(:);
    M(:,34)=y(:).*y(:).*x(:).*z(:);
    M(:,35)=z(:).*z(:).*x(:).*y(:);    
end
       
% estimate polynomial coefficients using least squares:
% M*a = data(:)
% M'*M*a = M'*data(:)
% a = inv(M'*M) * M'*data(:)
%cond(M'*M)

% Add a very small factor to the diagonal to Regularize when W is sparse 

a=(M'*Wm*M+diag(1e-10*ones(nterm,1)))\(M'*Wm*data(:));
%M*a - data(:)

%% estimate misfit.
Wnrm=W(:)./sum(W(:)); % sum(Wnrm)=1
msft=mean(abs((M*a - data(:)).*Wnrm))/mean(abs(data(:)));
   

end

function [H] = gaussian3d(siz,std) % Gaussian filter

 siz   = (siz-1)/2;
 [x,y,z] = meshgrid(-siz:siz,-siz:siz,-siz:siz);
 arg   = -(x.*x + y.*y + z.*z)/(2*std*std);

 H     = exp(arg);
 H(H<eps*max(H(:))) = 0;

 sumh = sum(H(:));
 if sumh ~= 0,
   H  = H/sumh;
 end;

end


function outputfigs(Ur,Ui,oss,ons,fold,opt)

% Output some figures
% opt=1 : output each stack as its own figure
% opt=2 : limited figures
% opt=3 : No figures

if(opt==1) % All outputs in thier own figure
    %montagestack(MagIm,[],'clip');title(['MagIm :: '  fold])
    montagestack(oss,[],'clip');title(['Octahedral Shear Strain :: '  fold]);colorbar
    montagestack(ons,[],'clip');title(['Octahedral Normal Strain :: '  fold]);colorbar

    montagestack(Ur(:,:,:,1));colormap(jet);title('Re x dsp')
    montagestack(Ui(:,:,:,1));colormap(jet);title('Im x dsp')
    montagestack(Ur(:,:,:,2));colormap(jet);title('Re y dsp')
    montagestack(Ui(:,:,:,2));colormap(jet);title('Im y dsp')
    montagestack(Ur(:,:,:,3));colormap(jet);title('Re z dsp')
    montagestack(Ui(:,:,:,3));colormap(jet);title('Im z dsp')
elseif(opt==2)
    montagestack(oss,[],'clip');title(['Octahedral Shear Strain :: '  fold]);colorbar    
end

end
    
    
    
    
    
    



