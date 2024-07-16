function [outfilename,Gp,Gdp] = run_LDI_fun_Smith(fitkernel,Delta,smoothfilter,nanfract,subjID,curldsp,motion,fms,voxelsize)
% runs LDI on a range of fitting kernels
% written by R. Okamoto 11/15/16 based on script


%   3) select parameters
%      a) fitting region (start w/ 3 x 3 x 3)
%      b) TLS parameter (0=ordinary least squares, 1000 fit x=f(y)
%   4) run TLS fitting of G' (Gp) and G'' (Gdp)
%   5) compute statistics (means, std, goodness of fit)

% inputs
% usecurl: 1 if using curl data, 0 if using displacement data
% fitkernel: vector of value for fitting,e.g. 3:2:9 (all values must be odd)
% Delta: noise parameter (lower value assumes variance in laplacian is large)  
% smoothfilter:  vector for smoothing displacement or curl data, 
% sd followed by kernel size e.g. [0.65 3];
% 
% 
% save data to the right file
if nargin < 5,
% LDIdir = uigetdir('.','Directory with curl or displacementinputs for LDI')
% cd(LDIdir)
end
% load unchanging fit data
%load ../allimg.mat subjID
% mydate is used if input files exist
% 
% if nargin < 6,
%     inputfile=uigetfile('*.mat','Select input file to use:');
% end
% if ~exist('mydate','var'), 
%     mydate=input('Enter date for fitting input data or in single quotes or <cr> for today:  ')
%     if isempty(mydate)
%         mydate=date;
%     end
% end
curlx = curldsp(:,:,:,1);
curly = curldsp(:,:,:,2);
curlz = curldsp(:,:,:,3);
% voxelsize = 2;
% fms = 50;
PC2micron = 2.6250;
rho = 1000;
h = voxelsize/1000;

todaydate=date;
% load(inputfile);
usecurl = 1;
if usecurl ==1,
%    eval(strrep(strrep('load #-$-curl_LDI_fitting_inputs.mat;','#',subjID),'$',mydate));
    savefiletempl=strrep(strrep('!_LDIresults_allslices-@-fit#_Delta$_smooth%_use.mat','!',subjID),'@',todaydate);
    curldata=cat(4,curlx,curly,curlz);
    nanmaskfinal=mkn(abs(curldata(:,:,:,1))>0);
elseif usecurl  == 0,
%    eval(strrep(strrep('load #-$-disp_LDI_fitting_inputs.mat;','#',subjID),'$',mydate));
    savefiletempl=strrep(strrep('!_LDIresults_allslices-@-fit#_Delta$_smooth%_use.mat','!',subjID),'@',todaydate);
    curldata=repmat(zeros(size(u)),[1 1 1 3]);%zeros(size(Brfo));  % make empty curl matrix
    nanmaskfinal=mkn(abs(u)>0);
    
else
    error('invalid flag')
end




%load([subjID,'_ExVivoMask','.mat'],'ExVivoMask');

% set parameters
%if nargin < 4, 
% Delta=1e-4
% nanfract=0.25

disp(['Delta:  ',num2str(Delta)])
disp(['nanfract:  ',num2str(nanfract)]) 
smoothL=0;
dvx=voxelsize;dvy=voxelsize;dvz=voxelsize;
fitrange= (fitkernel-1)/2;
%fitrange=1:7,
% display frequency and voxsize
disp(['Fitting at ',int2str(fms),' Hz w voxel sizes: ',num2str(dvx),', ',num2str(dvy),', ',num2str(dvx),' mm'])


% smoothfilter =[2 5];  %sd =1, window = 3;
% smoothfilter =[1 5];  % for 50 Hz
% smoothfilter=[];
% for Nxy=1:length(fitrange),
%     savefilename{Nxy}=strrep(savefiletempl,'#',int2str(2*fitrange(Nxy)+1));
%     savefilename{Nxy}=strrep(savefilename{Nxy},'$',strrep(sprintf('%6.4f',Delta),'.','p'));
%     if exist('smoothfilter'),
%         if ~isempty(smoothfilter)
%             savefilename{Nxy}=strrep(savefilename{Nxy},'%',[int2str(smoothfilter(2)),'sd',strrep(num2str(smoothfilter(1)),'.','p')]);
%         end
%     else
%         savefilename{Nxy}=strrep(savefilename{Nxy},'%',int2str(0));
%     end
%     if usecurl==0,
%         savefilename{Nxy}=strrep(savefilename{Nxy},'use','usedisp');
%     else
%         savefilename{Nxy}=strrep(savefilename{Nxy},'use','usecurl');
%     end
% end
%matlabpool(length(fitrange));
startTime=tic;
%parfor Nxy=fitrange,%  Nxy =1; % fit window is 2*Nxy + 1
u = motion(:,:,:,1);
v = motion(:,:,:,2);
w = motion(:,:,:,3);
for Nxy=1:length(fitrange)
    [Gp, Gdp,C, nre, R, dnre,smoothedD,smoothedL] = fitsmooth_uvw_rjofun_spiral_R2_parfor(cat(4,real(u),real(v),real(w)),cat(4,imag(u),imag(v),imag(w)),nanmaskfinal,voxelsize*[1 1 1], PC2micron*1e-3, fms, rho, h,Delta,nanfract,fitrange(Nxy), usecurl,curldata,smoothfilter,[]);
    % [Gp, Gdp,C, nre, R, dnre,smoothedD,smoothedL] = fitsmooth_uvw_rjofun_spiral_R2_parfor(cat(4,real(u),real(v),real(w)),cat(4,imag(u),imag(v),imag(w)),nanmaskfinal,voxelsize*[1 1 1], PC2micron*1e-3, fms, rho, h,Delta,nanfract,fitrange(Nxy), usecurl,curldata,smoothfilter,[],savefilename{Nxy});
end %parfor
%matlabpool close

elapsedTime=toc(startTime)

%%  summarize results
fitkernel=2*(fitrange)+1;
nanstat= @(x) [nanmean(x) nanstd(x) nanmedian(x) nanmin(x) nanmax(x)];
%summaryftempl=strrep(savefilename{end},['fit' int2str(fitkernel(end))],'#');
 for k=1:length(fitrange);
%eval(strrep('load $','$',strrep(summaryftempl,'#',['fit' int2str(fitkernel(k))])));
%load(savefilename{k})
Gpnan=Gp;Gpnan(Gp==0)=NaN;
Gdpnan=Gdp;Gdpnan(Gdp==0)=NaN;
nrenan=nre;nrenan(nre==0)=NaN;
Gpnan_all(:,:,:,k)=Gpnan; Gdpnan_all(:,:,:,k)=Gdpnan; nrenan_all(:,:,:,k)=nrenan;
LDIressum(k,:)=[nansum(~isnan(Gpnan(:))) nanstat(1e-3*Gpnan(:).*nanmaskfinal(:)) nanstat(1e-3*Gdpnan(:).*nanmaskfinal(:)) nanstat(nrenan(:))];
%LDIressum_ExVivoMask(k,:)=[nansum(~isnan(Gpnan(:).*mkn(ExVivoMask(:)))) nanstat(1e-3*Gpnan(:).*nanmaskfinal(:).*mkn(ExVivoMask(:))) nanstat(1e-3*Gdpnan(:).*nanmaskfinal(:).*mkn(ExVivoMask(:))) nanstat(nrenan(:).*mkn(ExVivoMask(:)))];
 end

outfilename=strrep(strrep(strrep('#_LDI_res_usecurl_%_$_summary.mat','#',subjID),'$',todaydate),'%',['smooth',int2str(smoothfilter(2)),'sd',strrep(num2str(smoothfilter(1)),'.','p')])
%save(outfilename,'Delta','UIUCmap','LDIressum', 'smoothfilter','smoothL','fitrange','nanmaskfinal','savefilename','Gpnan_all','Gdpnan_all','nrenan_all');
%save(outfilename,'LDIressum_ExVivoMask','-append');
 
