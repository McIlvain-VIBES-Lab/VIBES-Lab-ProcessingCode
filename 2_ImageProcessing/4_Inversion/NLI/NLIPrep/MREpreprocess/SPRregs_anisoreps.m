load registrations

combo_opt=1; 
s=size(subcortBINARY);
if(exist('subcortBINARY','var'))
    regsSC=subcortBINARY(:,:,:,1);
    for ii=2:size(subcortBINARY,4)
        r=subcortBINARY(:,:,:,ii);
        regsSC(r==1)=ii;
    end
end

if(exist('tractsBINARY','var'))
    regstracts=tractsBINARY(:,:,:,1);    
    for ii=2:size(tractsBINARY,4)
        r=tractsBINARY(:,:,:,ii);
        regstracts(r==1)=ii;
    end
end

if(exist('tractsLABELS','var'))
    regslabels=tractsLABELS(:,:,:,1);    
    for ii=1:size(tractsLABELS,4)
        r=tractsLABELS(:,:,:,ii);
        regslabels(r==1)=ii;
    end
end


if(combo_opt==1) % All Subcorts, some labels (including CC, fornix, CR, CST)
    rsclist=[1 2 3 4 5 6];
    rlablist=[3 4 5 6 7 8 23 24 25 26 27 28];
    rtractlist=[];
    
    
end

nreg=length(rsclist)+length(rlablist)+length(rtractlist);
regs=zeros([s(1) s(2) s(3)]);
ireg=0;
for ii=1:length(rsclist)
    ireg=ireg+1;
    r=subcortBINARY(:,:,:,ii);
    regs(r==1)=ireg;
end
for ii=1:length(rlablist)
    ireg=ireg+1;
    r=tractsLABELS(:,:,:,ii);
    regs(r==1)=ireg;
end
for ii=1:length(rsclist)
    ireg=ireg+1;
    r=tractsBINARY(:,:,:,ii);
    regs(r==1)=ireg;
end

save SPRregs.mat regs combo_opt rsclist rlablist rtractlist

