ls *.nod
nodf=input('node file name ','s');

junk=load(nodf);
x=junk(:,2);
y=junk(:,3);
z=junk(:,4);


mridim=[10 10 10];
xin=linspace(min(x)-0.05*range(x),max(x)+0.05*range(x),mridim(1));
yin=linspace(min(y)-0.05*range(y),max(y)+0.05*range(y),mridim(2));
zin=linspace(min(z)-0.05*range(z),max(z)+0.05*range(z),mridim(3));
regs=zeros(mridim);
regoutf=[nodf(1:end-4) '.regstack'];


fid=fopen(regoutf,'w');
fprintf(fid,'Array Dimensions \n'); 
fprintf(fid,'%7i %7i %7i \n',mridim);
fprintf(fid,'x coordinates \n');
fclose(fid);
dlmwrite(regoutf,xin,'-append','delimiter',' '); 
fid=fopen(regoutf,'a'); fprintf(fid,'y coordniates \n'); fclose(fid);
dlmwrite(regoutf,yin,'-append','delimiter',' ');
fid=fopen(regoutf,'a'); fprintf(fid,'z coordniates \n'); fclose(fid);
dlmwrite(regoutf,zin,'-append','delimiter',' ');
for ii=1:mridim(3)
    dlmwrite(regoutf,regs(:,:,ii),'-append','delimiter',' ');
end
