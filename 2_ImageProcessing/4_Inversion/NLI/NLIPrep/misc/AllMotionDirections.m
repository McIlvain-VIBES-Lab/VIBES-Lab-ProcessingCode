% Creates folders with all 8 unique combinations of directions in DirIndex

if(exist('HeaderData.mat','file'))
    load HeaderData.mat
    load MRE_3DMotionData.mat
    ftyp=1;
    Motion2Image=ones(3,3);
elseif(exist('NLI_MRE_Input.mat','file'))
    load NLI_MRE_Input
    ftyp=2;
end
load Mask.mat

fstm='test_alldirs_';

ii=0;

ii=ii+1;
d(ii).tag='100-010-001'
d(ii).dirs=[1 0 0;0 1 0;0 0 1]*Motion2Image;

ii=ii+1;
d(ii).tag='n100-010-001'
d(ii).dirs=[-1 0 0;0 1 0;0 0 1]*Motion2Image;

ii=ii+1;
d(ii).tag='100-0n10-001'
d(ii).dirs=[1 0 0;0 -1 0;0 0 1]*Motion2Image;

ii=ii+1;
d(ii).tag='n100-0n10-001'
d(ii).dirs=[-1 0 0;0 -1 0;0 0 1]*Motion2Image;

ii=ii+1;
d(ii).tag='010-100-001'
d(ii).dirs=[0 1 0;1 0 0;0 0 1]*Motion2Image;

ii=ii+1;
d(ii).tag='0n10-100-001'
d(ii).dirs=[0 -1 0;1 0 0;0 0 1]*Motion2Image;

ii=ii+1;
d(ii).tag='010-n100-001'
d(ii).dirs=[0 1 0;-1 0 0;0 0 1]*Motion2Image;

ii=ii+1;
d(ii).tag='0n10-n100-001'
d(ii).dirs=[0 -1 0;-1 0 0;0 0 1]*Motion2Image;


ndir=length(d);
cdir=pwd;
for ii=1:ndir
    mkdir([fstm d(ii).tag])
    if(ftyp==1)
        copyfile('MRE_3DMotionData.mat',[fstm d(ii).tag])
        DirIndex(1:3,4:6)=d(ii).dirs;
        save([fstm d(ii).tag '/HeaderData.mat'],'DirIndex','freqHz')
    elseif(ftyp==2)
        Motion2Image=d(ii).dirs;
        save ([fstm d(ii).tag '/NLI_MRE_Input.mat'],'Ur','Ui','Motion2Image','voxsize_mm','freqHz','AnatomicalMRI','RHcoord')
    end

        

    copyfile('Mask.mat',[fstm d(ii).tag])
end



