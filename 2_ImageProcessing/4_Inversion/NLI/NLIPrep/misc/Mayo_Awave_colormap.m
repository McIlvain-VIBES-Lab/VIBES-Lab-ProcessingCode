% Make awave colormap Mayo uses for displacements

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
colormap(gca,awave);
c=get(gca,'clim');
set(gca,'clim',[-max(abs(c)) max(abs(c))]);
