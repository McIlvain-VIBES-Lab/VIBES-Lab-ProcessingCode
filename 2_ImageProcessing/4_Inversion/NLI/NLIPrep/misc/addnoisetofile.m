% Function to add random gaussian noise to supplied rows of a file and
% output a seperate file
% Example addnoisetofile(dispfile.dsp,2:4,5) will create a file called
%         dispfile.dsp_5pcntnse with 5% noise added to columns 2 to 4
% Noise level based on the mean absolute value of all rows

function [newfname]=addnoisetofile(fname,cols,lev)

vals=load(fname);

newvals=vals;
mnv=mean(mean(abs(vals(:,cols))));
for jj=cols    
    newvals(:,jj)=vals(:,jj) + lev/100*mnv*randn(size(vals,1),1);
end

newfname=[fname '_' num2str(lev) 'pcntnse'];


% Create format string
% Check if first row is integer:
if(sum(mod(newvals(:,1),1))==0)
    fstr='%8i';
else
    fstr='%16.9e';
end

%add %16.9e for each row
for ii=2:size(vals,2)
    fstr=[fstr ' %16.9e'];
end


fid=fopen(newfname,'w');
fprintf(fid,[fstr ' \n'],newvals');
fclose(fid);

end
    
