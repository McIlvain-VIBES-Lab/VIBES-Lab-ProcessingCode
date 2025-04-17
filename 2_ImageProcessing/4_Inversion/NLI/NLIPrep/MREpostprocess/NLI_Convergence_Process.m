
% Store all iterations of all properties in a mat file

% filetime is the creation time in minutes. 

d=dir('*RE..*0001*01.mtr');

disp('Iteration 0001 Files in folder ::')
for ii=1:length(d)
    disp([int2str(ii) ' :: ' d(ii).name]);
end

flist=[];%1;
if(isempty(flist))
    flist=input('Datasets to process eg <[1 3 4]> ');
end

%QN3shell2p05isox_mtr1_YYZX_0DTInse_0dspnse.RE..0009.prop.05.mtr
%QN1shell2p05isox_mtr1_YYZX_0DTInse_0dspnse.reconind  
%QN4shell2p05isox_mtr1_YYZX_0DTInse_0dspnse.mtrmesh.03.nod 
for ii=1:length(flist)
    clear reconind nod convg filetime regconv globconv
    fname=d(flist(ii)).name;    
    
    % find the number of iterations
    nitr=1;
    fileexist=true;
    invfile=fname;
    while(fileexist)
        invfile(end-15:end-12)=sprintf('%4.4i',nitr);
        if(exist(invfile,'file'))
            nitr=nitr+1;
        else
            fileexist=false;
            nitr=nitr-1;
        end
    end
    
    % find the number of properties 
    numprop=1;
    propexist=true;
    invfile=fname;
    while(propexist)
        invfile(end-5:end-4)=sprintf('%2.2i',numprop);
        if(exist(invfile,'file'))
            numprop=numprop+1;
        else
            propexist=false;
            numprop=numprop-1;
        end
    end    
    disp(['Set : ' fname ' :: ' int2str(numprop) ' Properties ,' int2str(nitr) ' iterations'])
    %open the reconind file
    fid=fopen([fname(1:end-20) 'reconind']);
    reconind=fgetl(fid);
    reconind(reconind==' ')=[];
    reconind=reconind=='T';
    
    % load the nodefiles
    inod=1;
    fexist=true;
    globmaxreg=0;
    while(fexist)
        nodf=[fname(1:end-20) 'mtrmesh.' sprintf('%2.2i',inod) '.nod'];
        fexist=exist(nodf,'file');
        if(fexist)
            junk=load(nodf);
            nod(inod).xyz=junk(:,2:4);
            nod(inod).sense=junk(:,5);
            nod(inod).reg=junk(:,6);
            nod(inod).maxreg=max(nod(inod).reg);
            if(nod(inod).maxreg>globmaxreg)
                globmaxreg=nod(inod).maxreg;
            end
            nod(inod).ind=zeros(size(nod(inod).xyz));
            for kk=1:3
%                 nod(inod).res(kk)=mean(diff(unique(nod(inod).xyz(:,kk))));
%                 nod(inod).ind(:,kk)=round(nod(inod).xyz(:,kk)./nod(inod).res(kk));
%                 nod(inod).ind(:,kk)=nod(inod).ind(:,kk)-min(nod(inod).ind(:,kk))+2;
%                 nod(inod).siz(kk)=max(nod(inod).ind(:,kk))+1; % Size of array to give a 1 pixel buffer everywhere;            
%                 
                nod(inod).res(kk)=mean(diff(unique(nod(inod).xyz(:,kk))));
                nod(inod).ind(:,kk)=round((nod(inod).xyz(:,kk)-min(nod(inod).xyz(:,kk)))./nod(inod).res(kk));
                nod(inod).ind(:,kk)=nod(inod).ind(:,kk)+2;
                nod(inod).siz(kk)=max(nod(inod).ind(:,kk))+1; % Size of array to give a 1 pixel buffer everywhere;            
                
            end
            inod=inod+1;
        else
            numnodf=inod-1;
        end
    end
    
    % Preallocate the arrays
    filetime=zeros(nitr,1);
    regconv=zeros(nitr,numprop,2,globmaxreg,6); %Record the mean, stdev, max, min, 95th and 5th percentile for each region     
    globconv=zeros(nitr,numprop,2,6); %Record the mean, stdev, max, min, 95th and 5th percentile for the whole mesh;
    for itr=1:nitr
        for jprop=1:numprop
           % load the meshind file
           invfile=fname;
           invfile(end-5:end-4)=sprintf('%2.2i',jprop);
           invfile(end-15:end-12)=sprintf('%4.4i',itr);
           % Real part
           if(reconind(2*jprop-1))
               meshind=load([invfile(1:end-3) 'meshind']);
               convg(itr,jprop,1).stack=nan(nod(meshind).siz);
               convg(itr,jprop,1).mesh=meshind;               
           else
               convg(itr,jprop,1).stack=[];
           end
           % Imag part
           invfile(end-19:end-18)='IM';
           if(reconind(2*jprop))
               meshind=load([invfile(1:end-3) 'meshind']);
               convg(itr,jprop,2).stack=nan(nod(meshind).siz);
               convg(itr,jprop,2).mesh=meshind;
           else
               convg(itr,jprop,2).stack=[];
           end
        end
    end
    disp('Stacks preallocated')
    for itr=1:nitr
        disp(['Iteration ' int2str(itr) ' processing'])
        for jprop=1:numprop
           invfile=fname;
           invfile(end-5:end-4)=sprintf('%2.2i',jprop);
           invfile(end-15:end-12)=sprintf('%4.4i',itr);
           
           % Record the file creation time to gauge runtime
           if(jprop==1)
               % Accuracy of mins
%                 [dum,str] = unix(['ls -l ' invfile]);
%                 c=textscan(str,'%s');
%                 month=c{1}{6};
%                 day=c{1}{7};
%                 time=c{1}{8};
%                 formatin='dd-mmm-yyyy HH:MM';
%                 filetime(itr)=datenum([day '-' month '-2020 ' time],formatin)*24*60
                %accouracy of seconds
                [dum,str] = unix(['stat ' invfile]);
                c=textscan(str,'%s');
                day=c{1}{30}; % c{1}{29} should be 'modified'
                time=c{1}{31};
                formatin='yyyy-mm-dd HH:MM:SS';
                filetime(itr)=datenum([day ' ' time],formatin)*25*60;
           end
           
           % Real property
           if(reconind(2*jprop-1))                         
               % Load this property                   
               prop=load(invfile);
               prop=prop(:,2);     
               msh=convg(itr,jprop,1).mesh;
               convg(itr,jprop,1).stack(sub2ind(nod(msh).siz,nod(msh).ind(:,1),nod(msh).ind(:,2),nod(msh).ind(:,3)))=prop;  
               globconv(itr,jprop,1,1)=mean(prop(nod(msh).sense==1));
               globconv(itr,jprop,1,2)=std(prop(nod(msh).sense==1));
               globconv(itr,jprop,1,3)=min(prop(nod(msh).sense==1));
               globconv(itr,jprop,1,4)=max(prop(nod(msh).sense==1));
               globconv(itr,jprop,1,5)=prctile(prop(nod(msh).sense==1),5);
               globconv(itr,jprop,1,6)=prctile(prop(nod(msh).sense==1),95);
               
               for ireg=1:nod(msh).maxreg
                   if(~isempty(nod(msh).reg==ireg))
                       regconv(itr,jprop,1,ireg,1)=mean(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,1,ireg,2)=std(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,1,ireg,3)=min(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,1,ireg,4)=max(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,1,ireg,5)=prctile(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)),5);
                       regconv(itr,jprop,1,ireg,6)=prctile(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)),95);
                   else
                       regconv(itr,jprop,1,ireg,:)=NaN;
                   end
               end                      
           end
           % Imag property
           invfile(end-19:end-18)='IM';
           if(reconind(2*jprop))
               % Load this property
               prop=load(invfile);
               prop=prop(:,2);                       
               msh=convg(itr,jprop,2).mesh;
               convg(itr,jprop,2).stack(sub2ind(nod(msh).siz,nod(msh).ind(:,1),nod(msh).ind(:,2),nod(msh).ind(:,3)))=prop;   
               globconv(itr,jprop,2,1)=mean(prop(nod(msh).sense==1));
               globconv(itr,jprop,2,2)=std(prop(nod(msh).sense==1));
               globconv(itr,jprop,2,3)=min(prop(nod(msh).sense==1));
               globconv(itr,jprop,2,4)=max(prop(nod(msh).sense==1));
               globconv(itr,jprop,2,5)=prctile(prop(nod(msh).sense==1),5);
               globconv(itr,jprop,2,6)=prctile(prop(nod(msh).sense==1),95);
               
               for ireg=1:nod(msh).maxreg
                   if(~isempty(nod(msh).reg==ireg))
                       regconv(itr,jprop,2,ireg,1)=mean(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,2,ireg,2)=std(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,2,ireg,3)=min(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,2,ireg,4)=max(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)));
                       regconv(itr,jprop,2,ireg,5)=prctile(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)),5);
                       regconv(itr,jprop,2,ireg,6)=prctile(prop((nod(msh).sense==1)&(nod(msh).reg==ireg)),95);
                   else
                       regconv(itr,jprop,2,ireg,:)=NaN;
                   end
               end
           end          
        end
    end
    save([fname(1:end-21) '_' int2str(nitr) 'itr_convergence_stack.mat'],'reconind','nod','convg','filetime','regconv','globconv')
end




        
