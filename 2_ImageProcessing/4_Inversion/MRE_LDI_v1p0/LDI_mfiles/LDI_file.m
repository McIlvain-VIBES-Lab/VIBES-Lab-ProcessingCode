
% load ep2d_slch.mat
% ny = 50;
% nx = 50;
% nz = 20;
% ni = 3;
% Motion = cat(4,Xmotion,Ymotion,Zmotion);

% tspc = (2*pi/4);
% t= (0:tspc:((2*pi)-tspc))';
% tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
% ttmp = real(repmat(Motion,[1 1 1 1 4]).*exp(1i*tt));
% 
% PCdata = permute(ttmp,[1 2 3 5 4]);
tt = 0:(2*pi/8):(2*pi*(7/8))
Wdata = cat(5,Xmotion,Ymotion,Zmotion);
for ii = 1:length(tt)
PCdata(:,:,:,ii,:) = real(Wdata.*exp(-1i*tt(ii)));
end

userow = 1:size(PCdata,1);
usecol = 1:size(PCdata,2);
subjID = 'peyton_simulation_test';
nanmask = mask;
nanmask(nanmask == 0) = NaN;
nanmask = nanmask;
voxelsize = 1.5;
MAG = repmat(t2stack,[1 1 1 4 3]);
PC2micron = 2.6250;
fms = 60;
[Gp,Gdp] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID,fms);

Gc = Gp+1i*Gdp;
Mu = 2*(abs(Gc).^2)./(Gp+abs(Gc));
save('Mu.mat','Mu','Gp','Gdp')


