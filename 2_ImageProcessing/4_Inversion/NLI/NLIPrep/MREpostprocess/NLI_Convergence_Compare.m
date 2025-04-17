function NLI_Convergence_Compare(flist,truth,outfile,framerate,compopt,figpos)
% flist.name contains the convergence file from NLI_Convergence_Process
% flist.title (Optional) contains the titles for the plots
% flist.leg (optional) contains the legend titles. Defaults to first 3 characters of flist.name
% truthRE/truthIM contains the true values size(nx ny nz nprop 2) e.g. cat(5,propstack_RE, propstack_IM)
% plotopt=1 shows montages. Currently the only option.
% compopt=0 compares by iteration
% compopt>0 takes steps of compopt minutes and shows the closest image
% Example commands:
% flist=[dir('CG2*') dir('QN2*_YY*')];
% ii=1;flist(ii).title='CG2';flist(ii).leg=flist(ii).title;
% ii=2;flist(ii).title='QN2';flist(ii).leg=flist(ii).title;
% load('Shell_50Hzmod6_mtr1_xfacexdir.forward.props.mat')
% NLI_Convergence_Compare(flist,cat(5,propstack_RE,propstack_IM),'ConvergenceCG2vsQN2.avi',5)

plotopt=1; % for future options other than montagestack
clist='rgbcmykgbcmykgbcmykgbcmykgbcmykgbcmyk'; % Colors for each convergence trace
dotlist='.xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'; % Iteration markers for each region in regconv

if(nargin>2)
    % Ensure outfile has .avi at the end
    if(length(outfile)>4)
        if(~strcmp(outfile(end-3:end),'.avi'))
            outfile=[outfile '.avi'];
        end
    else
        outfile=[outfile '.avi'];
    end
else
    outfile='Convergence.avi';
end

if(nargin<4)
    framerate=5;
end

if(nargin<5)
    compopt=0;
end

if(nargin<6)
    figpos=[0.01 0.01 0.99 0.99];
end

v = VideoWriter(outfile);
v.FrameRate = framerate;
disp(['Output video file: ' outfile])
nf=length(flist);

maxitr=0;
maxtime=0;
for jj=1:length(flist)
    load(flist(jj).name);
    C(jj).convg=convg;
    C(jj).reconind=reconind;
    C(jj).filetime=filetime;
    C(jj).nrecon=sum(C(jj).reconind);
    C(jj).nitr=size(C(jj).convg,1);
    C(jj).nprop=size(C(jj).convg,2);
    C(jj).regconv=regconv;
    C(jj).globconv=globconv;
    
    if(C(jj).nitr>maxitr)
        maxitr=C(jj).nitr;
    end
    C(jj).filetime=C(jj).filetime-min(C(jj).filetime);
    if(C(jj).filetime(end)>maxtime)
        maxtime=C(jj).filetime(end);
    end
end

lastitr=zeros([size(C(1).convg(end,1,1).stack) C(1).nprop 2]);
for jprop=1:C(1).nprop
    for rr=1:2
        if(C(1).reconind(2*(jprop-1)+rr))
            lastitr(:,:,:,jprop,rr)=C(1).convg(end,jprop,rr).stack;
        end
    end
end
% if truth doesnt exist, take the last iteration of the first convergence
% set. Needs to be all the same mesh or this will crash. 
if((nargin<2)||isempty(truth))
    truth=lastitr;
    truthflag=false;
else
    truthflag=true;
end

% Define colorbar limits
cmax=zeros(size(truth,4),size(truth,5));
cmin=zeros(size(truth,4),size(truth,5));
for jprop=1:size(truth,4)
    for rr=1:2
        cmax(jprop,rr)=max(max(max(truth(:,:,:,jprop,rr))));
        cmin(jprop,rr)=min(min(min(truth(:,:,:,jprop,rr))));
    end
end

% cmax(1,1)=3800;
 cmax(1,2)=1500;



open(v);

curitr=zeros(nf,1);
if(compopt==0)
    numsteps=maxitr;
else
    numsteps=ceil(maxtime/compopt);
        
end
hf=figure('units','normalized','outerposition',figpos);
for istep=1:numsteps
    if(compopt~=0)
        tstep=istep*compopt;
    end
    for ifile=1:nf
        propcount=0;
        previtr=curitr;
        if(compopt~=0)
            delt=C(ifile).filetime-tstep;
            [junk,curitr(ifile)]=max(delt(delt<0));
        else
            curitr(ifile)=min(istep,C(ifile).nitr);
        end

        for iprop=1:C(ifile).nprop
            for rr=1:2
                if C(ifile).reconind(2*(iprop-1)+rr)
                    propcount=propcount+1;
                    if((istep==1)&&(ifile==1))
                        subplot(nf+2,C(ifile).nrecon,propcount);
                        montagestack(truth(:,:,:,iprop,rr),[6 11],'y','n');
                        caxis([cmin(iprop,rr) cmax(iprop,rr)]);
                        colorbar
                        if(propcount==1)
                            title('Real \mu')
                            ylabel('truth')
                        elseif(propcount==2)
                            title('Imag \mu')
                        elseif(propcount==3)
                            title('\phi')
                        elseif(propcount==4)
                            title('\zeta')
                        end
                            
                    end
                    
                    if((curitr(ifile)<=C(ifile).nitr)&&(curitr(ifile)~=previtr(ifile))) % DOnt update if past final iteration, or if iteration hasnt chaaaaaaaaaged with timestep
                        subplot(nf+2,C(ifile).nrecon,(ifile)*C(ifile).nrecon+propcount);
                        montagestack(C(ifile).convg(curitr(ifile),iprop,rr).stack,[],'n','n')
                        colorbar
                        caxis([cmin(iprop,rr) cmax(iprop,rr)]);
                    end
                                          
                    
                    if(propcount==1)
                        if(isfield(flist,'title')&&(~isempty(flist(ifile).title)))
                            ylabel(flist(ifile).title)                            
                        else
                            ylabel(flist(ifile).name(1:3))
                        end                        
                    end
                    
                    if(propcount==C(ifile).nrecon)
                        if(compopt==0)
                            ylabel([sprintf('%.0f',C(ifile).filetime(curitr(ifile))) ' mins'])                            
                        end
                    end
                    
                    subplot(nf+2,C(ifile).nrecon,(nf+1)*C(ifile).nrecon+propcount);
                    if(ifile==1)
                        hold off
                    end
                    for ireg=1:size(regconv,4)
                        h(ifile,propcount)=errorbar(1:C(ifile).nitr,C(ifile).regconv(:,iprop,rr,ireg,1),C(ifile).regconv(:,iprop,rr,ireg,2),[clist(ifile) dotlist(ireg)]);
                        hold on
                        plot(min(C(ifile).nitr,istep),squeeze(C(ifile).regconv(min(C(ifile).nitr,istep),iprop,rr,ireg,1)),'ko','linewidth',1.5) % Indicate current iteration
                    end
                    
                    
%                     plot(regconv(:,iprop,rr,:,3),'r-')
%                     plot(regconv(:,iprop,rr,:,4),'r-')
%                     plot(regconv(:,iprop,rr,:,5),'r:')
%                     plot(regconv(:,iprop,rr,:,6),'b:')
                    
                    
                    
                end
            end
        end
        
    end
    
    if(compopt==0)
        suptitle(['Iteration ' int2str(istep)])
    else
        suptitle(['Time ' int2str(tstep)])
    end
    %subplot(nf+2,C(ifile).nrecon,(nf+1)*C(ifile).nrecon);
    legend(h(:,1),flist(:).leg)
    
    pause(0.5)
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);









