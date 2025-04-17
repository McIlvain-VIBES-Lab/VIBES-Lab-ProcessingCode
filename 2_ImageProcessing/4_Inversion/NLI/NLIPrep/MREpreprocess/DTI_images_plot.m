function []=DTI_images_plot(DTIf,slc)
load(DTIf,'V1','FA')
if(nargin==1) % Plot montage
    figure;montage(permute(V1,[1 2 4 3]));
    title('Directions')
    
    FArep=zeros([size(FA,1) size(FA,2) 3 size(FA,3)]);
    for ii=1:3
        for jj=1:size(FA,3)
            FArep(:,:,ii,jj)=FA(:,:,jj);
        end
    end
    figure;montage(FArep.*permute(V1,[1 2 4 3]));
    title('colorFA')
    
elseif(nargin==2) % Plot slice
    figure;image(squeeze(V1(:,:,slc,:)));
    title('Directions')
    
    figure;image(repmat(FA(:,:,slc),[1 1 3]).*squeeze(V1(:,:,slc,:)));
    title('ColorFA')
end
    



