% Clear up MREv7inv diretory, leave only the latest iteration and some
% optional other iterations defined in svitr.
% Inputs: svitr = extra iteration to save. leave empty to save only the
% maximum iteration number.
% force = force delete, i.e. dont question. force='foce' will skip the
% questioning and just do it.
% Examples: CleanupMREv7invdirectory will save the latest iteration, and
% prompt before doing any deleteing
% CleanupMREv7invdirectory(10) will save iteration 10 and the latest iteration, and
% prompt before doing any deleteing
% CleanupMREv7invdirectory(10,'force') will save iteration 10 and the latest iteration, and
% offer no prompt before doing any deleteing
% CleanupMREv7invdirectory([],'force') will save the latest iteration, and
% offer no prompt before doing any deleteing
%
% Note that it will clean up ALL recon datasets in the dircctory
function CleanupMREv7p3invdirectory(svitr,force)

if(nargin==0)
    svitr=[]'
    force='no'
elseif(nargin==1)
    force='no'
end



% d1=dir('inv/*.dsp');
% d2=dir('inv/*.pre');


d3=dir('inv/*RE..0001.prop.01*mtr');

if(isempty(d3))
    error('No RE..0001 property files present, directory is already tidy')
end


% identify unique filestems
fstm(1).name=d3(1).name(1:end-21);
nstm=1;
for ii=2:length(d3)
    stm=d3(ii).name(1:end-21);
    unq=true;
    jj=0;
    while (unq&&jj<nstm)
        jj=jj+1;
        if(strcmp(fstm(jj).name,stm))
            unq=false;
        end
    end
    if(unq)
        nstm=nstm+1;
        fstm(nstm).name=stm;
    end
end

if(nargin==0)
    svitr=[];
    force='n';
end
if(nargin==1)
    force='n';
end

for istm=1:nstm
    pstm=dir(['inv/' fstm(istm).name '.RE..*.prop.01.mtr']);
    
    % Find maximum iteration number
    fstm(istm).mxitr=0;
    for ii=1:length(pstm)
        if(exist(['inv/' fstm(istm).name '.RE..' sprintf('%4.4i',ii) '.prop.01.mtr'],'file'))
            fstm(istm).mxitr=ii;
        end
    end
end
    
    
disp('oo')
disp(['oo ' int2str(length(fstm)) ' reconstructed datasets found'])
disp('oo')
for ii=1:length(fstm)
    disp([int2str(ii) ' :: ' fstm(ii).name ', Max iteration ' int2str(fstm(ii).mxitr)])
end

if (~strcmp(force,'force'))
    qn = input(['Remove all inversion files except ' int2str([svitr]) 'and maximum iteration, y(default)/n >>'],'s');
    if(isempty(qn)) 
        qn='y';
    end
else
    qn='y';
    disp(['Removing all inversion files except ' int2str([svitr mxitr]) '...'])
end

for ifle=1:length(fstm)
    if(strcmp(qn,'y'))
        for ii=1:fstm(ifle).mxitr-1
            if(sum([svitr fstm(ifle).mxitr]==ii)==0)
                %disp(ii)
                %ls(['inv/' fstm(ifle).name '.RE..' sprintf('%4.4i',ii) '.prop.*.mtr'])
                
                eval(['!rm inv/' fstm(ifle).name '.RE..' sprintf('%4.4i',ii) '.prop.*.mtr'])
                eval(['!rm inv/' fstm(ifle).name '.IM..' sprintf('%4.4i',ii) '.prop.*.mtr'])

                %if((~isempty(dir(['inv/' fstm(ifle).name '*' sprintf('%4.4i',ii) '*']))) && (~strcmp('mat',fstm(ifle).name(end-2:end))) &&  (~strcmp('fig',fstm(ifle).name(end-2:end))))
                    %eval(['!rm inv/' fstm(ifle).name '*.' sprintf('%4.4i',ii) '.*'])
                %else
                    %disp(['Keeping ' fstm(ifle).name])
                %end
            end
        end
    end
end
    
    