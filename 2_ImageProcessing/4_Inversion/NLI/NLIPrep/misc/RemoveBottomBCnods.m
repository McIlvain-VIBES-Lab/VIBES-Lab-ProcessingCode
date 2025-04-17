% Remove bottom boundary conditions from bnd.bcs file
% Run ~/code/misc/Bnod2bcs.x first to generate a file with all BCs. 
ls *bnd.bcs
bcf=input('BC file to remove bottom BCs from >> ','s');
ls *.nod
nodf=input('relevant node file >> ','s');


bc=load(bcf);
nod=load(nodf);
nod=nod(:,2:4);

bcnod=bc(:,2);
figure
plot3(nod(bcnod,1),nod(bcnod,2),nod(bcnod,3),'b.');
hold on

minz=min(nod(bcnod,3))+0.0001%0.0065;%0.000001;

btmbcind=nod(bcnod,3)<minz;
btmbcnods=bcnod(btmbcind);

plot3(nod(btmbcnods,1),nod(btmbcnods,2),nod(btmbcnods,3),'ro');
title('BC Nodes to be removed')

newbc=bc(~btmbcind,:);
newbc(:,1)=1:size(newbc,1);
newbcnod=newbc(:,2);
figure
plot3(nod(newbcnod,1),nod(newbcnod,2),nod(newbcnod,3),'g.')
title('Remaining BCs')
outfile=[bcf(1:end-3) 'nobtm.bcs'];

fid = fopen(outfile,'w');
fprintf(fid,'%7i  %7i  %i  %e  %e  %i  %e  %e  %i  %e  %e \n',newbc');
fclose(fid);
disp([int2str(size(bc,1)-size(newbc,1)) ' of ' int2str(size(bc,1)) ' BCs removed '])
disp(['Output file: ' outfile])