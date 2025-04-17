

load Rep01MRE_AP_2set_voxelmesh_G3300.v9p37NITI_isodamp_SPoffavlast08..0100.ReconProps.mat

Realmuref=Realmu;
MagImref=MagIm;
propstackref=cat(4,Realmu,Imagmu,Realphi,Realzeta,DR);

Realmureg=zeros([size(Realmu) 4]);
Imagmureg=zeros([size(Realmu) 4]);
Realphireg=zeros([size(Realmu) 4]);
Realzetareg=zeros([size(Realmu) 4]);
DRreg=zeros([size(Realmu) 4]);
MagImreg=zeros([size(Realmu) 4]);

Realmureg(:,:,:,1)=Realmu;
Imagmureg(:,:,:,1)=Imagmu;
Realphireg(:,:,:,1)=Realphi;
Realzetareg(:,:,:,1)=Realzeta;
DRreg(:,:,:,1)=DR;
MagImreg(:,:,:,1)=MagIm;




for ii=2:10
    load(['Rep' sprintf('%2.2i',ii) 'MRE_AP_2set_voxelmesh_G3300.v9p37NITI_isodamp_SPoffavlast08..0100.ReconProps.mat'])
    Realmui=Realmu;
    MagImi=MagIm;
    propstack=cat(4,Realmu,Imagmu,Realphi,Realzeta,DR);

    [propstack_reg,t2stack_reg]=Register_Stack_to_T2(MagImi,propstack,MagImref,voxsize_mm);
    
    
    Realmureg(:,:,:,ii)=Realmu;
    Imagmureg(:,:,:,ii)=Imagmu;
    Realphireg(:,:,:,ii)=Realphi;
    Realzetareg(:,:,:,ii)=Realzeta;
    DRreg(:,:,:,ii)=DR;
    MagImreg(:,:,:,ii)=MagIm;

    
    
end




%%
% Realmureg(isnan(Realmureg))=0;
% Imagmureg(isnan(Realmureg))=0;
% Realphireg(isnan(Realmureg))=0;
% Realzetareg(isnan(Realmureg))=0;
% DRreg(isnan(Realmureg))=0;
% MagImreg(isnan(Realmureg))=0;
close all
montagestack(mean(Realmureg/1000,4),[],'clip');title('Mean Realmu');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(std(Realmureg/1000,[],4),[],'clip');title('Std Realmu');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(mean(Imagmureg/1000,4),[],'clip');title('Mean Imagmu');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(std(Imagmureg/1000,[],4),[],'clip');title('Std Imagmu');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(mean(Realphireg,4),[],'clip');title('Mean Realphi');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(std(Realphireg,[],4),[],'clip');title('Std Realphi');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(mean(Realzetareg,4),[],'clip');title('Mean Realzeta');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(std(Realzetareg,[],4),[],'clip');title('Std Realzeta');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(mean(DRreg,4),[],'clip');title('Mean DR');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(std(DRreg,[],4),[],'clip');title('Std DR');colorbar;pause(0.1);drawnow;pause(0.1);

montagestack(mean(MagImreg,4),[],'clip');title('Mean T2');colorbar;pause(0.1);drawnow;pause(0.1);
montagestack(std(MagImreg,[],4),[],'clip');title('Std T2');colorbar;pause(0.1);drawnow;pause(0.1);


