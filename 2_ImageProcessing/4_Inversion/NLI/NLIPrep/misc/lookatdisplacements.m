%% Look at displacements


clear all
close all


%plot_typ=1; % Cycle through slices
%plot_typ=2; % all dsps 1 figure
plot_typ=3; % seperate figs ur ui vr vi wr wi
Mayo_Awave_colormap;
cmap=awave;
MREdef1='MRE_3DMotionData.mat';
MREdef2='NLI_MRE_Input.mat';
MREdef3='MRE_for_Inversion.mat';
maskf='Mask.mat'; % for ftyp 1 and 3

if(exist(MREdef1,'file'))
    matf=MREdef1;
    ftyp=1;
elseif(exist(MREdef2,'file'))
	matf=MREdef2;
    ftyp=2;
elseif(exist(MREdef3,'file'))
    matf=MREdef3;
    ftyp=3;
else
    ls *mat
    matf=input('Displacement mat file name >> ','s');
    ftyp=input('File type (1=MRE_3DMotionData, 2=NLI_MRE_Input, 3=mre_for_inversion) >> ');
end


load(matf)
if(ftyp<3)
    load(maskf)
end

if(ftyp==1)
    Ur=A.*cos(P);
    Ui=A.*sin(P);
    AnatomicalMRI=MagIm;
elseif(ftyp==3)
    Ur=zeros([size(Xmotion) 3]);
    Ui=zeros([size(Xmotion) 3]);
    Ur(:,:,:,1)=real(Xmotion);
    Ur(:,:,:,2)=real(Ymotion);
    Ur(:,:,:,3)=real(Zmotion);
    Ui(:,:,:,1)=imag(Xmotion);
    Ui(:,:,:,2)=imag(Ymotion);
    Ui(:,:,:,3)=imag(Zmotion);
    AnatomicalMRI=t2stack;
end
    


% Trim back A and P to the mask size
I=find(mask);
[i,j,k]=ind2sub(size(mask),I);
buf=2;
ilim=[max(min(i)-buf,1) : min(max(i)+buf,size(mask,1))];
jlim=[max(min(j)-buf,1) : min(max(j)+buf,size(mask,2))];
klim=[min(k) : max(k)];

mask=mask(ilim,jlim,klim);
AnatomicalMRI=AnatomicalMRI(ilim,jlim,klim);
Ur=Ur(ilim,jlim,klim,:);
Ui=Ui(ilim,jlim,klim,:);



for ii=1:3
    Ur(:,:,:,ii)=mask.*Ur(:,:,:,ii);
    Ui(:,:,:,ii)=mask.*Ui(:,:,:,ii);
end

titles(1,1:14)='U Displacement';
titles(2,1:14)='V Displacement';
titles(3,1:14)='W Displacement';

if(plot_typ==1)
for ii=1:size(Ur,3);
    
    h1=figure('units','normalized','position',[0,0,0.5,0.8]);
    subplot(2,2,1)
    imagesc(AnatomicalMRI(:,:,ii));axis equal; axis tight
    title('Anatomical')
    %set(gcf,'position',[116   500   560   420])
    h2=figure('units','normalized','position',[0.5,0,0.5,0.8])
    subplot(2,2,1)
    imagesc(AnatomicalMRI(:,:,ii));axis equal; axis tight
    title('Anatomical')
    %set(gcf,'position',[704   500   560   420])
    
        
    for jj=1:3
        figure(h1)
        subplot(2,2,jj+1)
        imagesc(Ur(:,:,ii,jj));axis equal; axis tight
        title(['Real ' titles(jj,:)])
        colormap(gca,cmap)
        colorbar
        figure(h2)
        subplot(2,2,jj+1)
        imagesc(Ui(:,:,ii,jj));axis equal; axis tight
        title(['Imag ' titles(jj,:)]) 
        colormap(cga,cmap)
        colorbar
    end
    figure(1)
    suptitle(['Slice ' int2str(ii) ' Press Space for next slice'])
    figure(2)
    suptitle(['Slice ' int2str(ii) ' Press Space for next slice'])
    pause
end
elseif(plot_typ==2)
% Create montage of all slices
asz=size(Ur);
alldsps=ones((asz(1)+1)*asz(3)+1,(asz(2)+1)*6+1)*max([Ur(:)' Ui(:)']);

% insert all displacement components into image
for ii=1:asz(3)
    for jj=1:3
        alldsps((ii-1)*(asz(1)+1)+1:(ii)*(asz(1)+1)-1,(jj-1)*(asz(2)+1)+1:(jj)*(asz(2)+1)-1)=Ur(:,:,ii,jj);
        alldsps((ii-1)*(asz(1)+1)+1:(ii)*(asz(1)+1)-1,3*(asz(2)+1)+(jj-1)*(asz(2)+1)+1:3*(asz(2)+1)+(jj)*(asz(2)+1)-1)=Ui(:,:,ii,jj);
    end
end
fsiz=18; % Font size
hf=figure;
ha=imagesc(alldsps);
colormap(gca,cmap);
axis equal;axis tight;hc=colorbar;
set(gca,'XTickLabel','','YTickLabel','')
set(hc,'FontSize',fsiz)
dirstr=pwd;slash=find(dirstr=='/');
ht=title(['Displacements, from folder ' dirstr(slash(end)+1:end)]);
set(ht,'FontSize',fsiz,'Interpreter','none')
hx=xlabel('Order: Re(u) Re(v) Re(w) Im(u) Im(v) Im(w)');
set(hx,'FontSize',fsiz,'Interpreter','none')
hy=ylabel(['Slice number  :   ' int2str(asz(3)) '<--  --> 1']);
set(hy,'FontSize',fsiz,'Interpreter','none')
elseif(plot_typ==3)
    montagestack(Ur(:,:,:,1));title('Real comp 1')
    colorbar;colormap(gca,cmap);saveas(gcf,'RealComp1.jpg','jpg')
    montagestack(Ur(:,:,:,2));title('Real comp 2')
    colorbar;colormap(gca,cmap);saveas(gcf,'RealComp2.jpg','jpg')
    montagestack(Ur(:,:,:,3));title('Real comp 3')
    colorbar;colormap(gca,cmap);saveas(gcf,'RealComp3.jpg','jpg')
    montagestack(Ui(:,:,:,1));title('Imag comp 1')
    colorbar;colormap(gca,cmap);saveas(gcf,'ImagComp1.jpg','jpg')
    montagestack(Ui(:,:,:,2));title('Imag comp 2')
    colorbar;colormap(gca,cmap);saveas(gcf,'ImagComp2.jpg','jpg')
    montagestack(Ui(:,:,:,3));title('Imag comp 3')
    colorbar;colormap(gca,cmap);saveas(gcf,'ImagComp3.jpg','jpg')
end
    





    