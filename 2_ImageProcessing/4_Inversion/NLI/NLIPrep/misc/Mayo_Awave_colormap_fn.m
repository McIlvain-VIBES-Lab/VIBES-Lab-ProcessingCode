% Make awave colormap Mayo uses for displacements
function Mayo_Awave_colormap_fn(H,clim)
% Changes figure to mayo awave colormap
% Inputs: H = axis handle
%         clim(optional) value to use for max and min of colormap
% eg Mayo_Awave_colormap_fn(gca,3) will set the current axis to the mayo
%                                  colormap with clim = [-3 3]

zto1=0:0.01:1;
onetoz=1:-0.01:0;
one=ones(size(zto1));
z=zeros(size(zto1));

%    Cyan  ->  Blue  -> Black  ->  Red  -> Yellow
% R   0         0         0         1         1
% G   1         0         0         0         1
% B   1         1         0         0         0

awaver=[z z zto1 one];
awaveg=[onetoz z z zto1];
awaveb=[one onetoz z z ];

awave=[awaver' awaveg' awaveb'];
% Use colormap(awave) to use this colormap.

% Apply this colormap to gca

colormap(H,awave);

if(nargin==1)
    c=get(gca,'clim');    
elseif(nargin==2)
    if(length(clim)==1)
        c=[-clim clim];
    elseif(length(clim==2))
        c=clim;
    end
else
    error(['nargin = ' int2str(nargin)])
end

set(gca,'clim',[-max(abs(c)) max(abs(c))]);