% Plot a montage a stack of images, autoscaled with optional size argument.
% Inputs: stackin = 3d image stack
%         monsize = optional montage dimensions. Set to [] for default
%         clip = (optional) set to 'y' or 'clip' image to only plot 
%           non-zero data from the stack
%         newfig = (optional). set to 'n' to plot in the currently active figure
%           window, otherwise a new figure will be opened.          
% Examples: montagestack(A(:,:,:,2)) will plot the y directed displacement
%             amplitudes 
%           montagestack(structname(3).vals,[4 3],'y') will plot the 3rd stack
%             in structname, in a 4x3 montage, with the images clipped to 
%             only show non-zero data.
%           montagestack(A,[],'y') will plot the clipped montage of A,
%             using the default automatically calculated montage size.
%           montagestack(mask,[6 2]) will display a 6x2 montage of the
%             mask, without clipping.
%           montagestack(A,[],'y','n') will display a clipped montage of A,
%             and will not create a new figure window
function [h]=montagestack(stackin,monsiz,clip,newfig)

stackin(isnan(stackin))=0;

if(nargin==1)
    monsiz=[];
end

if(nargin<3)
    clip='n';
end

if(nargin<4)
  newfig='y';
end

if (strcmp(clip,'y')||strcmp(clip,'clip'))
    buf=2; % Number of pixels for buffer
    jnk=find(stackin);
    [v1,v2,v3]=ind2sub(size(stackin),jnk);
    l1mn=max(min(v1)-buf,1);
    l1mx=min(max(v1)+buf,size(stackin,1));
    l2mn=max(min(v2)-buf,1);
    l2mx=min(max(v2)+buf,size(stackin,2));
else
    l1mn=1;
    l1mx=size(stackin,1);
    l2mn=1;
    l2mx=size(stackin,2);
end
stack=stackin(l1mn:l1mx,l2mn:l2mx,:);

dim=size(stack);

if(~strcmp(newfig,'n'))
  figure
end

if(nargin==1||isempty(monsiz))
    [h]=montage(reshape(stack,[dim(1) dim(2) 1 dim(3)]),'displayrange',[]);    
elseif(nargin>1)
    [h]=montage(reshape(stack,[dim(1) dim(2) 1 dim(3)]),'displayrange',[],'Size',monsiz);
else
    error(['Incorrect number of arguments :: nargin = ' int2str(nargin)])
end
drawnow
end
