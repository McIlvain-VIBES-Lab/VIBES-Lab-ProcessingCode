% MRE_process_folder: Processes whole folder of MRE reconstructions.

% Assumes only 1 mesh inside 'hex' folder
% Runs the convergence postprocessor if required
% Processes all reconstructions in the inv folder, and interpolates the
% first 2 mtr files (assuming that is the final iteration and averaged
% file. 





%close all;clear all;clc

function MRE_v7p3_process_folder(fold)


if(nargin<1) % Default is current folder 
  fold=[];  
elseif(fold(end)~='/')
    fold(end+1)='/';    
end

if exist('hex','dir')&&exist('tet','dir')
    meshtyp=input('hex and tet folders present, which one to process (default hex) >> ','s');
    if(~strcmp(meshtyp,'tet'))
        meshtyp='hex';
    end
elseif exist('hex','dir')
    meshtyp='hex';
elseif exist('tet','dir')
    meshtyp='tet';
end

cd([fold meshtyp])
f=dir;
cd([f(end).name '/inv'])
pwd
%% Find mtr.0001 to postprocess

f0001=dir('*RE..0001.prop.01.mtr');
for ii=1:length(f0001)
    disp(['Convergence postprocessing : ' f0001(ii).name])
    eval(['!~/Discovery_home/code/convcode/MREpostprocessv7p3.x ' f0001(ii).name])
      
end

fconv=dir(['*.prop.01.RE.convinfo']);
for ii=1:length(fconv)
    Plotpostprocess7p3(fconv(ii).name,'fig')  
end

fnod=dir('*01.nod');

for ii=1:length(fnod)
    disp(['Interpolation : ' fnod(ii).name])
    nodstm=fnod(ii).name(1:end-14);
    fmtr=dir([nodstm(1:end-1) '*RE.*prop.01.mtr']);
    for jj=1:2
        disp(['Mtr file ' fmtr(jj).name])
        MRE_plotv7p3_mask_Linear_Function(fnod(ii).name,fmtr(jj).name)
    end
end

!cp *.fig *.mat ~/Discovery_home/data/matfiletransfer
if(isempty(fold))
    cd('../../../')
else
    cd('../../../..')
end

end

% if(nargin<1)
%     d=dir([nodstm(1:end-1) '*RE.*prop.01.mtr']);
%     if length(d)~=0
%         mfind=false;
%         for ii=1:length(d)
%             disp(['File ' int2str(ii) ' :: ' d(ii).name])
%         end
%     else
%         d=dir([nodstm(1:end-1) '*RE.*prop.01.mf.mtr']);
%         mfind=true;
%         for ii=1:length(d)
%             disp(['File ' int2str(ii) ' :: ' d(ii).name])
%         end
%     end
%     clear n
%     n=input('Number of file to use (default = last)  >> ');
%     if(isempty(n))
%         n=length(d);
%     end
%     
%     if ~mfind
%         mtrstm=d(n).name(1:end-20);
%         itrstr=d(n).name(end-15:end-12);
%     else
%         mtrstm=d(n).name(1:end-23);
%         itrstr=d(n).name(end-18:end-15);
%     end
%     
% else
%     if(strcmp(mtrf(end-6:end),'mf.mtr'))
%         mfind=true;
%     else
%         mfind=false;
%     end
%     
%     if ~mfind
%         mtrstm=mtrf(1:end-20);
%         itrstr=mtrf(end-15:end-12);
%     else
%         mtrstm=mtrf(1:end-23);
%         itrstr=mtrf(end-18:end-15);
%     end    
% end    
    


