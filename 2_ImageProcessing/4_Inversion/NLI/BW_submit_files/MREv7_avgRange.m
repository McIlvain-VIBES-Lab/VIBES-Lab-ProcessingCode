function MREv7_avgRange(N_avstart,N_avend)
% MRE-v7-Convergence

% Plots Some convergence information For MRE-v7
%clear all
%close all
mf=false;
[f]=getfilenamemfDef('qn',['*.RE.*0001.prop.01'],mf);
if ~mf
    pref=f(1:end-20);
    mid='.prop.';
else
    pref=f(1:end-23);
    mid='.prop.';
end

% Determine number of properties
for ii=1:99
    ps=sprintf('%2.2i',ii);
    fn=[pref 'RE..' '0001' mid ps '.mtr'];
    if(exist(fn,'file'))
        np=ii;
        suff{ii}='.mtr';
        mfpropind(ii)=false;
    else
        fn=[pref 'RE..' '0001' mid ps '.mf.mtr'];
        if(exist(fn,'file'))
            np=ii;
            suff{ii}='.mf.mtr';
            mfpropind(ii)=true;
        end
    end
end

% Determine number of Iterations
for ii=1:999
    istr=sprintf('%4.4i',ii);
    fn=[pref 'RE..' istr mid '01' suff{1}];
    if(exist(fn,'file'))
        nitr=ii;
        if(ii==999)
            disp('WARNING:: currently only 999 iterations can be detected')
        end
    end    
end

disp(['Plotting ' int2str(np) ' Properties and ' int2str(nitr) ' Iterations.'])

% Load recon index
rcnindf=[pref 'reconind'];
if(exist(rcnindf,'file'))
    fid=fopen(rcnindf);
    fr=fgetl(fid);
    fr=fr(fr~=' ');
    rcnind=(fr=='T');
else % Use default rcnind
    rcnind=[true true false true true false];
end


%% Prompt to average values
%avstart=(nitr-Navg+1);
avstart=N_avstart;
avend=N_avend;
%
if(isempty(avstart))
   disp('No averaging performed')
else
    sstr=sprintf('%4.4i',avstart);
    for ip=1:np
        for jj=1:2
            if(jj==1)
                RI='RE..';
                RI_t='Real';
            elseif(jj==2)
                RI='IM..';
                RI_t='Imag';
            end
            if(rcnind(2*(ip-1)+jj)) 
                Nav=(avend-avstart+1);
                for ii=avstart:avend
                    istr=sprintf('%4.4i',ii);
                    ps=sprintf('%2.2i',ip);
                    fn=[pref RI istr mid ps suff{ip}];
                    junk=load(fn);                
                    if(ii==avstart)
                        if ~mfpropind(ip)
                            propval=junk(:,2);
                            SD=junk(:,2).^2;
                        else
                            propval=junk(:,2:3);
                            SD=junk(:,2:3).^2;
                        end
                    else
                        if ~mfpropind(ip)
                            propval=propval+junk(:,2);
                            SD=SD + junk(:,2).^2;
                        else
                            propval=propval+junk(:,2:3);
                            SD=SD + junk(:,2:3).^2;
                        end
                    end
                end
                propval=propval./Nav;
                SD=sqrt( Nav/(Nav-1)*(1/(Nav)*SD - propval.^2));
                % Write averaged file
                avfnme=[pref RI 'avfrom' sstr '.' istr mid ps suff{ip}];
                fid=fopen(avfnme,'w');
                if ~mfpropind(ip)
                    fprintf(fid,'%7i %12.5e\n',[(1:size(propval,1))' propval]');
                else
                    fprintf(fid,'%7i %12.5e %12.5e\n',[(1:size(propval,1))' propval]');
                end
                fclose(fid);
                disp(['Averaged property file ' avfnme ' created'])
                
                avfnme=[pref RI 'stdfrom' sstr '.' istr mid ps suff{ip}];
                fid=fopen(avfnme,'w');
                if ~mfpropind(ip)
                    fprintf(fid,'%7i %12.5e\n',[(1:size(propval,1))' SD]');
                else
                    fprintf(fid,'%7i %12.5e %12.5e\n',[(1:size(propval,1))' SD]');
                end
                fclose(fid);
                disp(['Property standard deviation file ' avfnme ' created'])
            else % Prop not reconstructed, copy final file
                istr=sprintf('%4.4i',avend);
                ps=sprintf('%2.2i',ip);
                fn=[pref RI istr mid ps suff{ip}];
                avfnme=[pref RI 'avfrom' sstr '.' istr mid ps suff{ip}];
                junk=load(fn);
                fid=fopen(avfnme,'w');
                if ~mfpropind(ip)
                    fprintf(fid,'%7i %12.5e\n',junk');
                else
                    fprintf(fid,'%7i %12.5e %12.5e\n',junk');
                end
                fclose(fid);
                disp(['Property file ' avfnme ' copied from final iteration'])
                
                avfnme=[pref RI 'stdfrom' sstr '.' istr mid ps suff{ip}];
                junk(:,2)=0; % Zero std for unreconstructed params
                fid=fopen(avfnme,'w');
                if ~mfpropind(ip)
                    fprintf(fid,'%7i %12.5e\n',junk');
                else
                    fprintf(fid,'%7i %12.5e %12.5e\n',junk');
                end
                fclose(fid);
                disp(['Property standard deviation file ' avfnme ' created'])
            end
        end
    end
end            
            
            
