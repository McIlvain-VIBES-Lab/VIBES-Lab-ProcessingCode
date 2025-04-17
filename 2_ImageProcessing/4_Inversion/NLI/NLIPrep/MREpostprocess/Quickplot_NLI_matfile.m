
% Quickly plot all the niti paramter images from a NLI mat file.
% Just the isotropic damping version

function Quickplot_NLI_matfile(matf)

if (nargin==0) % Prompt for file
    [matf]=getfilename('NLI .mat file to plot images ','*.mat')
end

load(matf);
mtyp=invparam.modeltype;

if(mtyp==1)
    msk=RealShear>0;
    I=find(msk);
    [i,j,k]=ind2sub(size(msk),I);
    klim=[min(k):max(k)];
    DR(isnan(DR))=0;
    montagestack(RealShear(:,:,klim)/1000,[],'clip');title('Real mu');colorbar;drawnow;pause(0.01)
    montagestack(ImagShear(:,:,klim)/1000,[],'clip');title('Imag mu');colorbar;drawnow;pause(0.01)
    montagestack(msk(:,:,klim).*DR(:,:,klim),[],'clip');title('Damping Ratio');colorbar;drawnow;pause(0.01)

elseif(mtyp==3) % Visco comp
    msk=RealShear>0;
    I=find(msk);
    [i,j,k]=ind2sub(size(msk),I);
    klim=[min(k):max(k)];
    DR(isnan(DR))=0;
    montagestack(RealShear(:,:,klim)/1000,[],'clip');title('Real mu');colorbar;drawnow;pause(0.01)
    montagestack(ImagShear(:,:,klim)/1000,[],'clip');title('Imag mu');colorbar;drawnow;pause(0.01)
    montagestack(msk(:,:,klim).*DR(:,:,klim),[],'clip');title('Damping Ratio');colorbar;drawnow;pause(0.01)
    montagestack(RealLambda(:,:,klim)/1000,[],'clip');title('Real lambda');colorbar;drawnow;pause(0.01)
    
elseif(mtyp==4) % Poro
    msk=RealShear>0;
    I=find(msk);
    [i,j,k]=ind2sub(size(msk),I);
    klim=[min(k):max(k)];
    montagestack(RealShear(:,:,klim)/1000,[],'clip');title('Real mu');colorbar;drawnow;pause(0.01)
    montagestack(RealLambda(:,:,klim)/1000,[],'clip');title('Real Lambda');colorbar;drawnow;pause(0.01)
    montagestack(RealHC(:,:,klim),[],'clip');title('log10(HC)');colorbar;drawnow;pause(0.01)
    
    
elseif(mtyp==6)
    msk=Realmu>0;
    montagestack(Realmu/1000,[],'clip');title('Real mu');colorbar;drawnow;pause(0.01)
    montagestack(Imagmu/1000,[],'clip');title('Imag mu');colorbar;drawnow;pause(0.01)
    DR=Imagmu./Realmu/2;
    DR(isnan(DR))=0;
    montagestack(DR,[],'clip');title('Damping Ratio');colorbar;drawnow;pause(0.01)
    montagestack(Realphi,[],'clip');title('phi');colorbar;drawnow;pause(0.01)
    montagestack(Realzeta,[],'clip');title('zeta');colorbar;drawnow;pause(0.01)
end

end

