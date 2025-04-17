function [maxsiz,vol]=AnalyzeMesh(nodf,elmf,maskf)
%ANALYZEMESH Summary of this function goes here
%   Detailed explanation goes here

nod=load(nodf);
x=nod(:,2);
y=nod(:,3);
z=nod(:,4);

nn=length(x);

elm=load(elmf);
if(size(elm,2)<=7) % Probably a tet mesh
    in=elm(:,2:5);
    npe=4;
elseif(size(elm,2)>27) % Probably a hex27 mesh
    in=elm(:,2:28);
    npe=27;
end
ne=size(elm,1);

maxsiz=zeros(ne,1);
vol=zeros(ne,1);
lengthratio=zeros(ne,1);
maxlength=zeros(ne,1);
minlength=zeros(ne,1);
for el=1:ne
    xe=x(in(el,:));
    ye=y(in(el,:));
    ze=z(in(el,:));
    
    if(npe==4)
        vol(el)=1/6*det([xe ye ze ones(4,1)]);
        edgelength(1)=sqrt((xe(1)-xe(2))^2+(ye(1)-ye(2))^2+(ze(1)-ze(2))^2);
        edgelength(2)=sqrt((xe(1)-xe(3))^2+(ye(1)-ye(3))^2+(ze(1)-ze(3))^2);
        edgelength(3)=sqrt((xe(1)-xe(4))^2+(ye(1)-ye(4))^2+(ze(1)-ze(4))^2);
        edgelength(4)=sqrt((xe(2)-xe(3))^2+(ye(2)-ye(3))^2+(ze(2)-ze(3))^2);
        edgelength(5)=sqrt((xe(2)-xe(4))^2+(ye(2)-ye(4))^2+(ze(2)-ze(4))^2);
        edgelength(6)=sqrt((xe(3)-xe(4))^2+(ye(3)-ye(4))^2+(ze(3)-ze(4))^2);
        maxlength(el)=max(edgelength);
        minlength(el)=min(edgelength);
        lengthratio(el)=minlength(el)/maxlength(el);        
    end    
    maxsiz(el)=max([range(xe) range(ye) range(ze)]);
end

if(nargin>=3)
    load(maskf)
    disp(['Mask used to create mesh has ' int2str(sum(mask(:))) ' points'])
end
disp(['Mesh has ' int2str(nn) ' nodes and ' int2str(ne) ' elements'])
disp(['Max edge length =' num2str(mean(maxlength)) ' +/- ' num2str(std(maxlength))]) 
disp(['Edge length ratio =' num2str(mean(lengthratio)) ' +/- ' num2str(std(lengthratio))]) 



% Plot some histograms
figure;histogram(vol);
title('Element volume distribution')
figure;histogram(maxsiz);
title('Element max size')
figure;histogram(lengthratio);
title('min edge/max edge')
figure;histogram(maxlength);
title('Element max length')
figure;histogram(minlength);
title('Element min length')



