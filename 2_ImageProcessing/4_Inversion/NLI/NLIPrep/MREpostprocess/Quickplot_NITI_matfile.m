
% Quickly plot all the niti paramter images from a NLI mat file.
% Just the isotropic damping version

function Quickplot_NITI_matfile(matf)

if (nargin==0) % Prompt for file
    [matf]=getfilename('NLI .mat file to plot images ','*.mat')
end

load(matf);
mtyp=invparam.modeltype;

if(mtyp==1)
    montagestack(RealShear/1000,[],'clip');title('Real mu');colorbar;drawnow;pause(0.01);saveas(gcf,'InvisoRealMu.jpg','jpg')
    montagestack(ImagShear/1000,[],'clip');title('Imag mu');colorbar;drawnow;pause(0.01);saveas(gcf,'InvisoImagMu.jpg','jpg')
    montagestack(DR,[],'clip');title('Damping Ratio');colorbar;drawnow;pause(0.01);saveas(gcf,'InvisoDR.jpg','jpg')
    

elseif(mtyp==6)
montagestack(Realmu/1000,[],'clip');title('Real mu');colorbar;drawnow;pause(0.01);saveas(gcf,'InvanisoRealMu.jpg','jpg')
montagestack(Imagmu/1000,[],'clip');title('Imag mu');colorbar;drawnow;pause(0.01);saveas(gcf,'InvanisoImagMu.jpg','jpg')

DR=Imagmu./Realmu/2;DR(isnan(DR))=0;
montagestack(DR,[],'clip');title('Damping Ratio');colorbar;drawnow;pause(0.01);saveas(gcf,'InvanisoDR','jpg')
montagestack(Realphi,[],'clip');title('phi');colorbar;drawnow;pause(0.01);saveas(gcf,'Invanisophi.jpg','jpg')
montagestack(Realzeta,[],'clip');title('zeta');colorbar;drawnow;pause(0.01);saveas(gcf,'Invanisozeta.jpg','jpg')
end

end

