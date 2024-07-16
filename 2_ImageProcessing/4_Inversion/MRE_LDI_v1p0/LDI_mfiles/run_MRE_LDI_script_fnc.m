function [Gp,Gdp] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID,fms)
%script to run LDI for spiral or EPI MRE
% written by R. Okamoto 
% updated 11/22/16
%
% This script assumes you have unwrapped displacement data (encoded as phase) stored in a 5D
% array where array dimension 4 is time (usually 8 time points0 and array
% dimension 5 is displacement component (x, y, z) where x corresponds to
% the direction of column of your data volume, y corresponds to the row and
% z corresponds to the slice (or stack) direction

%it also assumes you have made a mask of your data (probably needed for
%unwrapping. The size of the mask is the same as the first three dimensions
% of the phase contrast data
%this script is best run in sections (denoted by %% at the beginning)

%% Step 1 
% load data file with unwrapped phase contrast data and other variables
% clear all; %otherwise things will get messed up

% [MREdatafile,MREdir]=uigetfile('Select file containing MRE data:  ');
% cd(MREdir);
% load(MREdatafile);
if ~exist('PCdata','var') % then using our WUSTL naming convention
    PCdata=PCWU;clear PCWU;
    save(MREdatafile); % this removes PCWU from the MRE datafile
end
[nrow, ncol, nslice,ntp,ndir]=size(PCdata);
if ~exist('voxelsize','var')
    if exist('voxsize','var'),
        voxelsize=voxsize;
    else
       voxelsize=input('Enter isotropic voxel size in mm:  ');
    end
    save(MREdatafile,'voxelsize','-append');
end
%% Step 2:  Remove Bulk Motion

% BMRflag = input('Do you want to remove the bulk motion? 1=yes, 0 = no?  ');
BMRflag = 0;
if BMRflag==1,
    clear PCdata_BMR;
    for k=1:8
        [PCdata_BMR(:,:,:,:,k),~,X,Y,Z, coef_PCdata_BMR(:,k)]=calc_3Ddisp_BulkMotionRemoved(squeeze(PCdata(:,:,:,k,:)).*repmat(nanmask,[1 1 1 3]));
    end
    
    PCdata_BMR=permute(PCdata_BMR,[1 2 3 5 4]);
    
    
    %get bulk motion by PCdata-PCdata_BMR
    PCdata_BMR_nonan=PCdata_BMR;
    PCdata_BMR_nonan(isnan(PCdata_BMR_nonan))=0;
    
    
    PCdata1=PCdata_BMR_nonan;  %run OSS_SNR with the new PCdata_BMR
    figure();
    plot(coef_PCdata_BMR'); title('Bulk motion coefs');legend('rot z','-rot y','rot x','x trans','y trans','z trans')
    %save(MREdatafile,'PCdata_BMR','BMRflag','coef_PCdata_BMR', '-append')  % save unwrapping info
    
    
    % calculate bulk motion at center of masked volume
    [XX,YY,ZZ]=meshgrid(1:nrow,1:ncol,1:nslice);
    xcent=round(sum(mask(:).*XX(:))./sum(mask(:)));
    ycent=round(sum(mask(:).*YY(:))./sum(mask(:)));
    zcent=round(sum(mask(:).*ZZ(:))./sum(mask(:)));
    loc=[ycent xcent zcent];
    
    %  u_RB = u_R2*Y + u_R3*Z + u_R4
    %  v_RB = -u_R2*X +v_R3*Z + v_R4
    %  w_RB = -u_R3*X - v_R3*Y + w_R4
    % matrix is assembled so that the 6 variables can be found using least
    % squares
    %  solution vector is {u_R2 u_R3 vR3 u_R4 v_R4 w_R4}
    % coef 1 is rotation about Z
    % coef 2 is negative rotation about Y
    % coef 3 is rotation about X
    % units of coef 1 to 3 are radians/per VOXEL (not mm)
    % units of coef 4 to 6 are phase contrast radians
    coef_PCdata_BMR_loc=coef_PCdata_BMR;
    
    for i=1:8
        coef_PCdata_BMR_loc(4,i) = coef_PCdata_BMR(4,i)+coef_PCdata_BMR(1,i)*loc(2)+coef_PCdata_BMR(2,i)*loc(3);
        coef_PCdata_BMR_loc(5,i) = coef_PCdata_BMR(5,i)-coef_PCdata_BMR(1,i)*loc(1)+coef_PCdata_BMR(3,i)*loc(3);
        coef_PCdata_BMR_loc(6,i) = coef_PCdata_BMR(6,i)-coef_PCdata_BMR(2,i)*loc(1)-coef_PCdata_BMR(3,i)*loc(2);
    end
    save(MREdatafile,'coef_PCdata_BMR_loc','loc','-append');
    
    disp('Location of Bulk Motion Coefficients:')
    disp(loc);
    disp('Bulk Motion Coefficients at this location:')
    disp(coef_PCdata_BMR_loc);
    % NOTE: to get bulk motion at any voxel, compute PCdata_Bulk
    PCdata_Bulk=PCdata-PCdata_BMR;
end % if BMRflag
%% 
disp('Step 3: select and prepare data to use for fitting');

% BMRtype = input('Select PC data to use for fitting, enter 1 for PCdata, 2 for PCdata_BMR:  ');
BMRtype = 1;
if BMRtype ==1,
    P =PCdata;
elseif BMRtype ==2,
    P = PCdata_BMR;
else
    error('Could not find data matrix')
end
%plot data
% plotyes=input('Do you want to plot the data, 1 if yes, 0 if no:  ');
plotyes = 0;
if plotyes,
    vis5d_profile(P,colormap(jet),'PC data for fitting')
    input('Enter any key to continue :  ');
end
% remove any nans
if sum(isnan(P(:)))>0,
    P(isnan(P))=0;
end
if ~exist('subjID','var'),
    subjID=input('Enter subject ID in single quotes:  ');
end
%check for PC2micron
if ~exist('PC2micron','var')
    PC2micron=input('No PC2micron found in input file, enter value here:  ')
end
% set values for userow and usecol if not set already
if ~exist('userow','var'),
    if exist('cropdim','var'),
        tmp=strfind(cropdim,',');
        eval(['userow=',cropdim(1:(tmp(1)-1)),';']);
        eval(['usecol=',cropdim((tmp(1)+1):(tmp(2)-1)),';']);
    else
        userow=1:size(P,1);
        usecol=1:size(PCdata,2);
    end
end
% select slices to process
%useslice = input('Enter vector of slices to use for fitting, or CR to use all slices:  ');
useslice = [];
if isempty(useslice),
    useslice=1:nslice;
end
disp(['using slices ',int2str(useslice(1)),' to ',int2str(useslice(end))])
% define region to fit - important
P=P(userow,usecol,useslice,:,:); 
cropmask=mask(userow,usecol,useslice);
cropnanmask=nanmask(userow,usecol,useslice);
disp('Size of data volume to fit: ')
size(cropmask)
%% Step 4: calculate strain and filtered displacement
% parameters from Matt
filtwidth=0.5; % Width of Gaussian smoothing filter (pixels). Wider = more smoothing.     
ord=2; % order of polynomial to fit (only 2nd/quadratic 3rd/cubic or 4th/quartic currently supported.
Nfit=3; % Size of block to fit data to. Needs to be odd. Must be 6*x filtwidth or greater

% filtwidth=0.75; % to match our method
% ord=2; % order of polynomial to fit (only 2nd/quadratic 3rd/cubic or 4th/quartic currently supported.
% Nfit=5; % Size of block to fit data to. Needs to be odd. Must be 6*x filtwidth or greater
disp(['Filtering using filtwidth of : ',num2str(filtwidth),' pixels, poly order = ',int2str(ord),',filter size of ',int2str(Nfit), ' voxels']);
OSS_SNRfname=[subjID,'_OSS_SNR'];
if BMRtype==2,
    OSS_SNRfname=strrep(OSS_SNRfname,'_OSS','_BMR_OSS');
end
if length(useslice)<nslice,
    OSS_SNRfname=strrep(OSS_SNRfname,'_OSS',['_slices',int2str(useslice(1)),'to',int2str(useslice(end)),'_OSS']);
end
disp(['Saving data to: ',OSS_SNRfname])
% using McGarry code, outputs updated on 6/6/15

[OSS_SNR,Motion_SNR,OSS_SNR_Dist,Motion_SNR_Dist,oss,ons,curldsp,sall,Urf,Uif,MSNR_4D]=Strain_SNR_from_phase_MMG_rjoedit(P*PC2micron/1e6,cropmask,voxelsize*1e-3*[1 1 1],filtwidth,ord,Nfit,OSS_SNRfname);

% define a median function for nans 
stats= @(x) x(round(sum(~isnan(x))*[0.1 0.5 0.9]))';
mynanmed=@(x) nanmedian(x(:));
% find median strains
disp('Shear strains, then normal strains, then curl components, then OSS then ONS\n')
disp('10%, 50%, 90% of range\n')
uselen=length(userow)*length(usecol)*length(useslice);
for k=1:6,
    medstr(k)=mynanmed(abs(sall(:,:,:,k)).*mkn(abs(sall(:,:,:,k))>0))*1e6;
    strainstatsMGG(k,:)=stats(sort(reshape(abs(sall(:,:,:,k)).*mkn(abs(sall(:,:,:,k))>0),[uselen 1])))*1e6;
end
strainstatsMGG=strainstatsMGG([4:6 1:3],:);% reorder to match
for k=1:3,
    medcurl(k)=mynanmed(abs(curldsp(:,:,:,k)).*mkn(abs(curldsp(:,:,:,k))>0))*1e6;
    strainstatsMGG(k+6,:)=stats(sort(reshape(abs(curldsp(:,:,:,k)).*mkn(abs(curldsp(:,:,:,k))>0),[uselen 1])))*1e6;
end
strainstatsMGG(10,:)=stats(sort(reshape(oss.*mkn(oss>0),[uselen 1])))*1e6;
strainstatsMGG(11,:)=stats(sort(reshape(ons.*mkn(ons>0),[uselen 1])))*1e6

disp([' Number of datapoints:   ', int2str(sum(ons(:)>0))])
disp('OSS_SNR Motion_SNR')
[OSS_SNR;Motion_SNR]
[UU,EM]=Errormap_MMG_rjoedit(P);  % this calculates the error map and the harmonic of the displacements before ANY smoothing
RelUnc=2/size(P,4)*EM./abs(UU);  % this is to be consistent w/ Matt's Strain_SNR_from_phase code. Note RelUnc is the inverse of MotionSNR components
RMSRelUnc=sqrt(sum(RelUnc.^2,4));  % compute RMS of uncertainty of three motion components assuming they are independent
% RMS Rel Unc should measure the overall reliability of a data point. Since
% neighboring datapoints are used in the calculation of strain or curl, neighboring Rel Unc is also important
stmp=sort(OSS_SNR_Dist(:).*mkn(OSS_SNR_Dist(:)>0));
rtmp=sort(RMSRelUnc(:).*mkn(RMSRelUnc(:)>0));
% put these lines in cells below OSS_SNR and Motion_SNR
[stmp(round([0.1 0.5 0.9]*sum(~isnan(stmp))))';rtmp(round([0.1 0.5 0.9]*sum(~isnan(rtmp))))']



% resize Urf and Uif to remove padding
lf=floor(Nfit/2); %Limits either side of centre
bufsiz=lf+1; % Size of zero buffer
Urf=Urf((bufsiz+1):(end-bufsiz),(bufsiz+1):(end-bufsiz),(bufsiz+1):(end-bufsiz),:);
Uif=Uif((bufsiz+1):(end-bufsiz),(bufsiz+1):(end-bufsiz),(bufsiz+1):(end-bufsiz),:);

% check size
if any((size(squeeze(P(:,:,:,1,:)))-size(Urf))~=0),
    error('Inconsistent sizes');
end


%% Step 5
sprintf('step 5, plot data if desired \n\n')
% plotyes = input('Do you want to plot filtered displacement and curl? 1=yes, 0 = no:   ');
plotyes = 0;
if plotyes,
    if ~exist('Urf','var'),
        whos
        RealDispVar = input('Enter variable with real displacements: ');
        ImagDispVar = input('Enter variable with imag displacements: ');
    else
        RealDispVar='Urf';ImagDispVar='Uif';
    end
    eval(strrep('RealDispVar=$;','$',RealDispVar));
    eval(strrep('ImagDispVar=$;','$',ImagDispVar));

%     %find
%     if ~exist('userow','var'),
%         if exist('cropdim','var'),
%             tmp=strfind(cropdim,',');
%             eval(['usecol=',cropdim((tmp(1)+1):(tmp(2)-1)),';']);
%             eval(['userow=',cropdim(1:(tmp(1)-1)),';']);
%         else
%             [userow,usecol]=size(Urf)
%         end
%     end
    % plot data - this command uses a 5D viewing tool
    h1 =vis5d_profile(1e6*mkn(cat(5,RealDispVar ,ImagDispVar )),colormap(jet),'Re/Im disp in microns',[-1 1]*max(abs(1e6*RealDispVar(:))),{'u/v/w','Re/Im'});%% step 2
    %fitorder and kernel size are chosen to match Matt McGarry's values; strain values will differ due to greater smoothing
    %compared to Matt
   
    h2 = vis5d_profile(mkn(cat(5,real(curldsp),imag(curldsp))),colormap(jet),'Re/im curl in mm/mm',[-1 1]*max(abs(curldsp(:))),{'curl x/y/z','Re/Im'});
    go=input('Hit enter to continue: ');
end

%% Step 6
sprintf('\n step 6: setup fitting parameters \n')

% LDIdir = uigetdir('.','Directory to use for LDI')
% cd(LDIdir)
% LDIdir=pwd;

h = 1e-3*voxelsize;
% dens = input('Enter density in g/cm^3:  ');
dens = 1;
rho = 1e3*dens;  % density of porcine brain
UIUCmap=UIUCelastogram;
   
disp([ 'Voxel size = ',num2str(h),' in m'])
disp(['density = ',num2str(rho),'kg/m^3'])


%fms = input('Enter freq in Hz:  ');
%fms = 50;
disp(['Using freq of:  ',num2str(fms), '  Hz']) 
    

% smooth2=0; % do not smooth Laplacian
% umask = 0; % mask before smoothing?
% 
% % TLS settings for all fit size
% nanfract = 0.5;
% 
% usecurl=input('Enter 1 to use curl data, or 0 to use displacements:  ');
usecurl=1;
nanmask=mkn(mask);  % make a mask of nans

disp(['Using voxel size of: ',num2str(h*1e3), '  mm']) 

% define mm/voxel scaling in mm
    dvx=h*1e3;dvy=h*1e3;dvz=h*1e3;
    voxdims=[dvx dvy dvz];

% convert PC to mm instead of microns
    disp(['Using PC2micron of: ',num2str(PC2micron), '  microns/rad']) 
%     PC2mm=PC2micron*1e-3;
    
    
%
dispdir={'u','v','w'};
for k=1:3,
    eval(strrep('$=Urf(:,:,:,k)+1i*Uif(:,:,:,k);','$',dispdir{k}));
end
if usecurl==1,
    curlx = curldsp(:,:,:,1);
    curly = curldsp(:,:,:,2);
    curlz = curldsp(:,:,:,3);
end

% save the result
varlist={'fms','cropmask','cropnanmask','userow','usecol','useslice','PC2micron','u','v','w','subjID','voxelsize','usecurl','rho','UIUCmap','h'};
%cd(LDIdir);
if usecurl==1,
 fitinputfilename=[subjID,'-',date,'-','curl_LDI_fitting_inputs.mat'];
 varlist(length(varlist)+(1:3))={'curlx','curly','curlz'};
 save(fitinputfilename,varlist{:})
else
  fitinputfilename=[subjID,'-',date,'-','disp_LDI_fitting_inputs.mat'];
   save(fitinputfilename,varlist{:})
end


%% Step 7
sprintf('\n step 7: finalize variable parameters and run LDI \n');

% Delta=input('Enter delta value or <cr> for default value of 0.01: ');
% if isempty(Delta),
%     Delta=0.01;
% end
Delta = 0.1; % default 0.01; changed 180518

% fitkernel=input('Enter range for fitting kernel (vector of odd numbers) or <cr> for default value of 3:2:9 :  ');
% if isempty(fitkernel),
%     fitkernel=3:2:9;
% end

% fitkernel = 3;
% fitkernel = 5;
 fitkernel = 7;
% default values

% smoothfilter = input('Enter smoothing filter sd and kernel size as vector or <cr> for default [1 5]:  ');
% if isempty(smoothfilter),
% smoothfilter = [ 1 5];
% end

smoothfilter = [1 5];

disp(['Smoothing filter sd = ',int2str(smoothfilter(1)),',  kern = ',int2str(smoothfilter(2))])
nanfract =0.5;

%run_LDI_fun(fitkernel,Delta,smoothfilter,nanfract)

motion = cat(4,u,v,w);
[outfilename,Gp,Gdp] = run_LDI_fun_Smith(fitkernel,Delta,smoothfilter,nanfract,subjID,curldsp,motion,fms,voxelsize);

% LDI_ShearModulus_Plot(outfilename);

% load(outfilename,'LDIressum')
% disp('Number of datapoints: ')
% fprintf('%d\n',LDIressum(:,1));
% disp('Real(G) mean/std/median/min/max')
% fprintf('%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n',LDIressum(:,2:6)')
% disp('Imag(G) mean/std/median/min/max')
% fprintf('%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n',LDIressum(:,7:11)')
% disp('normalized residual error (nre) mean/std/median/min/max')
% fprintf('%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n',LDIressum(:,12:16)')

%clear all