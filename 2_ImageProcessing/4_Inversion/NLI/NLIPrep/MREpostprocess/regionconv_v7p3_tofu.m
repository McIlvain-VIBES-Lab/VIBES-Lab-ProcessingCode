
% Plot the convergence of each region with iteration
clear all;close all
%fstm='v7p32_SF0p008';


ii=0;

%% QN
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed5en9_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en9_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed5en10_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en10_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed5en11_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en11_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed5en12_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en12_nodelay_SFoff_v7p32_QN5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed5en13_nodelay_SFoff_v7p32_QN5';

%% CG
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en9_nodelay_SFoff_v7p32_CG5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en10_nodelay_SFoff_v7p32_CG5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en11_nodelay_SFoff_v7p32_CG5';
% ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en12_nodelay_SFoff_v7p32_CG5';

%ii=ii+1;fstm(ii).name='inv/TofuT18gel100_222_SPfixed1en13_nodelay_SFoff_v7p32_QN5';

% ii=ii+1;fstm(ii).name='QN_SP1en11_fixed';
% ii=ii+1;fstm(ii).name='QN_SP1en12_fixed';
% 
% ii=ii+1;fstm(ii).name='CG_SP5en9_fixed';
% ii=ii+1;fstm(ii).name='CG_SP1en9_fixed';
% ii=ii+1;fstm(ii).name='CG_SP5en10_fixed';
% ii=ii+1;fstm(ii).name='CG_SP1en10_fixed';
% ii=ii+1;fstm(ii).name='CG_SP5en11_fixed';
% ii=ii+1;fstm(ii).name='CG_SP1en11_fixed';
% ii=ii+1;fstm(ii).name='CG_SP1en12_fixed';

% Large Seg

% ii=ii+1;fstm(ii).name='inv/T18_100_SP1en9_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP5en10_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP1en10_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP5en11_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP1en11_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP5en12_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP1en12_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP5en13_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP1en13_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SP1en14_SFoff_largeseg_v7.32.visc';
% ii=ii+1;fstm(ii).name='inv/T18_100_SPoff_SF0.0015_largeseg_v7.32.visc';
% SP=[1e-9 5e-10 1e-10 5e-11 1e-11 5e-12 1e-12 5e-13 1e-13 1e-14 1e-15];
% Large Mis-Seg

% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP1en9_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP5en10_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP1en10_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP5en11_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP1en11_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP5en12_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP1en12_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP5en13_SF_off_QN5';
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP1en13_SF_off_QN5';
% 
% ii=ii+1;fstm(ii).name='inv/T18_100_large_misseg_SP1en14_SF_off_QN5';
% SP=[1e-9 5e-10 1e-10 5e-11 1e-11 5e-12 1e-12 5e-13 1e-13 1e-14];

%ii=ii+1;fstm(ii).name='COMET-101_SPR_voxelmesh_G3300.v7.3.inv.iso.incomp.visc_SPon';


d=dir('*.01.nod');
ii=1;fstm(1).name=d(1).name(1:end-15);
SP=[1e-9];
nsets=length(fstm);
disp(['Filename : ' fstm(1).name])

nprop=3;
cols=['bo';'rx';'g^';'cv';'ms';'k.'];


for nn=1:nsets

rcnf=[fstm(nn).name '.reconind'];
mshindf=[fstm(nn).name '.meshind'];
basisf=[fstm(nn).name '.basis'];


fid=fopen(rcnf);
f=fgetl(fid);
f=f(f~=' ');
rcnind=f=='T';


meshind=load(mshindf);
basis=load(basisf);

elmbasis=false;
if(max(basis(:))==2)
    elmbasis=true;
end


%% Read Mesh info
d=dir([fstm(nn).name '.mtrmesh*nod']);
nmesh=length(d);
maxreg=0;
for ii=1:nmesh
    junk=load(['inv/' d(ii).name]);
    mesh(ii).nod(:,:)=junk(:,2:4);
    
    if(size(junk,2)>4)
        mesh(ii).nodsense=junk(:,5);
    else
        mesh(ii).nodsense=ones(size(junk,1),1);
    end
    
    if(size(junk,2)>5)
        mesh(ii).nodreg=junk(:,6);
    else
        mesh(ii).nodreg=zeros(size(junk,1),1);        
    end
    maxreg=max(maxreg,max(mesh(ii).nodreg.*mesh(ii).nodsense));
    
    
    junk=load(['inv/' d(ii).name(1:end-3) 'elm']);
    mesh(ii).in=junk(:,2:9);
    mesh(ii).ne=length(junk);

    if(size(junk,2)>28)
        mesh(ii).elmsense=junk(:,29);
    else
        mesh(ii).elmsense=ones(size(junk,1),1);
    end

    if(size(junk,2)>29)
        mesh(ii).elmreg=junk(:,30);
    else
        mesh(ii).elmreg=zeros(size(junk,1),1);
    end
    maxreg=max(maxreg,max(mesh(ii).elmreg.*mesh(ii).elmsense));
    
    if(elmbasis)
        % Find element centroids to locate material property points
        for jj=1:mesh(ii).ne
            for kk=1:3
                mesh(ii).elmcent(jj,kk)=mean( mesh(ii).nod(mesh(ii).in(jj,:),kk) );
            end
        end
    end
    
end


%% Read Property files

mid='.prop.';
suff='.mtr';
% Determine number of Iterations
for ii=1:999
    istr=sprintf('%4.4i',ii);
    fn=[fstm(nn).name '.RE..' istr mid '01' suff];
    if(exist(fn,'file'))
        nitr=ii;
        if(ii==999)
            disp('WARNING:: currently only 999 iterations can be detected')
        end
    end    
end


meanvals=zeros([nprop,2,nitr,maxreg]);
stdvals=zeros([nprop,2,nitr,maxreg]);



for ii=1:nprop
    propstr=sprintf('%2.2i',ii);
    for jj=1:2
        if(jj==1)
            reim='.RE..';
        else
            reim='.IM..';
        end
        msh=meshind(jj,ii);
           
        
        if(rcnind(2*(ii-1)+jj))
            
            for kk=1:nitr
                istr=sprintf('%4.4i',kk);
                fn=[fstm(nn).name reim istr mid propstr suff];
                junk=load(fn);
                vals=junk(:,2);
                
                for rr=1:maxreg
                    if(basis(ii)==1)
                        reg=mesh(msh).nodreg.*mesh(msh).nodsense==rr;
                    else
                        reg=mesh(msh).nodreg.*mesh(msh).nodsense==rr;
                    end
                    
                    meanvals(ii,jj,kk,rr)=mean(vals(reg));
                    stdvals(ii,jj,kk,rr)=std(vals(reg));
                end
            end
            
            figure
            errorbar(1:nitr,meanvals(ii,jj,:,1),stdvals(ii,jj,:,1),cols(1,:))
            hold on
            for rr=2:maxreg
                errorbar(1:nitr,meanvals(ii,jj,:,rr),stdvals(ii,jj,:,rr),cols(rr,:))
            end
            
            title([fstm(nn).name 'Property ' int2str(ii) '-' reim(2:3)],'Interpreter','none')
            xlabel('iteration')
            ylabel('Property value')
            saveas(gcf,[fstm(nn).name 'prop' int2str(ii) reim(2:3) 'regconv.fig'],'fig')
            
        end
    end
end

save([fstm(nn).name '.regconv.mat'],'meanvals','nitr','stdvals','cols','reim','maxreg')

bvalsR(nn)=meanvals(1,1,end,1);
i1valsR(nn)=meanvals(1,1,end,2);
i2valsR(nn)=meanvals(1,1,end,3);

bvalsI(nn)=meanvals(1,2,end,1);
i1valsI(nn)=meanvals(1,2,end,2);
i2valsI(nn)=meanvals(1,2,end,3);
end
%%
%SP=[5e-9 1e-9 5e-10 1e-10 5e-11 1e-11 5e-12 1e-12 5e-13];
%SP=[1e-9 1e-10 1e-11 1e-12];

save('Regional_conv_rightseg.mat','SP','bvalsR','i1valsR','i2valsR','bvalsI','i1valsI','i2valsI')


figure
semilogx(SP,bvalsR(:),'bx-',SP,i1valsR(:),'rx-',SP,i2valsR(:),'gx-')
legend('QN Background','QN lower inclusion','QN upper Inclusion')
title('Real Shear Mod Final Iteration')

figure
semilogx(SP,bvalsI(:),'bx-',SP,i1valsI(:),'rx-',SP,i2valsI(:),'gx-')
legend('QN Background','QN lower inclusion','QN upper inclusion')
title('Imag Shear Mod Final Iteration')
