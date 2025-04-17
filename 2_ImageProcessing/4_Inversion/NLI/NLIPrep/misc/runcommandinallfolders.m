% Run a command in all folders
% Useful Examples
% Plot the shear modulus in each folder and save as a tif image
% runcommandinallfolders('load tetresults;load simulationinfo;montagestack(stacks(1).vals);colormap(jet);caxis([0 13000]);colorbar;title([''freq = '' int2str(freqHz)]);saveas(gcf,[''LE_3p5_freq'' int2str(freqHz) ''.tif'']) ')
% clean up all the excess iterations in the folders, saving the 77th
% iteration
% runcommandinallfolders('CleanupINVpardirectory(77)')
% Can split up into multiple commands, anything with an ! mark wont work
% with other matlab commands
% runcommandinallfolders('cd MRPE;','!qsub MRPE-submit','cd ..')
function runcommandinallfolders(varargin)
% Changes into each directory and runs the commands specified as inputs 'cmd'
% eg runcommandinallfolders('!qsub MPICH2par-submit') will run 

ninps=size(varargin,2);
dlist=dir;



nf=length(dlist);

for ii=3:nf
    if(isdir(dlist(ii).name));
        cd(dlist(ii).name);
        for jj=1:ninps
            cmd=varargin{jj};
            disp(['Running :: ' cmd ' in directory '  dlist(ii).name])
            eval(cmd);
        end
        cd ..
    end
end