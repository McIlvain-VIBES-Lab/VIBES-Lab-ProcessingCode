function moco3D_meanshot

% read_shift = 0;
% phase_shift = -35/240;
% % phase_shift = 0;

load recon_info/corners.mat

load recon_info/kspace_info.mat
load ascconv.mat

nsl = ascconv.sKSpace.lPartitions;
ntp = ascconv.lRepetitions +1;
nc = size(ascconv.asCoilSelectMeas.asList,2);
nl = ascconv.sKSpace.lPhaseEncodingLines;
nslab = ascconv.sGroupArray.asGroup.nSize;


%% Perform motion correction

fprintf(1,'Performing motion correction\n')

nav_vol_tmp = zeros(nroLR, LR_SIZE_Z, nl,nc,nsl,nslab,ntp);
for kk = 1:ntp
    for ii = 1:nslab
        load(sprintf('raw_data/raw_nav_%d_%d',ii,kk))
        nav_vol_tmp(:,:,:,:,:,ii,kk) = nav_vol;
    end
end
clear nav_vol

nav_imgb0 = zeros(LR_SIZE_XY, LR_SIZE_XY, LR_SIZE_Z, nc);
nav_imgb = zeros(LR_SIZE_XY, LR_SIZE_XY, LR_SIZE_Z, nc);
t_nav0=zeros(LR_SIZE_XY,LR_SIZE_XY,LR_SIZE_Z);
t_nav=zeros(LR_SIZE_XY,LR_SIZE_XY,LR_SIZE_Z);
nav_mean = zeros(nroLR, LR_SIZE_Z,nc);
mask = ones(LR_SIZE_XY, LR_SIZE_XY, LR_SIZE_Z);

dkx0 = zeros(nc, 1);
dky0 = zeros(nc, 1);
dkz0 = zeros( nc, 1);
dPhi0 = zeros(nc, 1);

dkx = zeros(nc, 1);
dky = zeros(nc, 1);
dkz = zeros(nc, 1);
dPhi = zeros(nc, 1);
corrections = zeros(4,nl, nsl, nslab, ntp);


for kk = 2:ntp
    for ii = 1:nslab
        % FIND BEST FIT SHOT
        for cc = 1:nc
            nshot = nl*nsl;
            nav_shots = reshape(nav_vol_tmp(:,:,:,cc,:,ii,kk),nroLR,LR_SIZE_Z,nshot);
            nav_mean(:,:,cc) = mean(nav_shots(:,:,(2:(nshot-1))),3);
            
%             for ss = 1:nshot
%                 nav_err(ss,cc) = norm(col(nav_shots(:,:,ss)-nav_mean),2);
%             end
        end
        
%         nav_err_tot = sum(nav_err,2);
%         nav_ind = find(nav_err_tot == min(nav_err_tot),1);
%         
%         [nl_ind,nsl_ind] = ind2sub([nl nsl],nav_ind);
        
        
        for cc = 1:nc
            % REFERENCE
            % cent_k = find(abs(nav_vol_tmp(:, :, nl_ind,cc,nsl_ind,ii,kk)) == max(abs(col(nav_vol_tmp(:, :,nl_ind,cc,nsl_ind,ii,kk)))), 1);
            cent_k = find(abs(nav_mean(:, :, cc)) == max(abs(col(nav_mean(:, :,cc)))), 1);
            [i_xy i_z] = ind2sub(size(nav_mean(:, :, cc)), cent_k);
            % store the shifts and offsets
            dkx0(cc) = kxLR(i_xy);
            dky0(cc) = kyLR(i_xy);
            dkz0(cc) = kzLR(i_z);
            dPhi0(cc) = angle(nav_mean(i_xy,i_z,cc).*exp(-1i*(read_shift.*kxLR(i_xy)+phase_shift.*kyLR(i_xy)+2/(2*LR_SIZE_Z)*kzLR(i_z))*2*pi));
            
            for idxn = 1:LR_SIZE_Z
                t_nav0(:, :, idxn) = k2image(kxLR, kyLR, nav_mean(:, idxn, cc).*exp(-1i*(read_shift.*kxLR+phase_shift.*kyLR+2/(2*LR_SIZE_Z)*kzLR(idxn))*2*pi), wwLR, LR_SIZE_XY, 2.5);
            end
            nav_imgb0(:,:,:,cc) = fftshift(ifft(ifftshift(t_nav0, 3), LR_SIZE_Z, 3),3);
        end
        
        
        for jj = 1:nsl
            for ll = 1:nl
                
                for cc = 1:nc
                    % NON-REFERENCE
                    cent_k = find(abs(nav_vol_tmp(:, :, ll,cc,jj,ii,kk)) == max(abs(col(nav_vol_tmp(:, :,ll, cc,jj,ii,kk)))), 1);
                    [i_xy i_z] = ind2sub(size(nav_vol_tmp(:, :,ll, cc,jj,ii,kk)), cent_k);
                    % store the shifts and offsets
                    dkx(cc) = kxLR(i_xy);
                    dky(cc) = kyLR(i_xy);
                    dkz(cc) = kzLR(i_z);
                    dPhi(cc) = angle(nav_vol_tmp(i_xy,i_z, ll, cc,jj,ii,kk).*exp(-1i*(read_shift.*kxLR(i_xy)+phase_shift.*kyLR(i_xy)+2/(2*LR_SIZE_Z)*kzLR(i_z))*2*pi));
                    
                    for idxn = 1:LR_SIZE_Z
                        t_nav(:, :, idxn) = k2image(kxLR, kyLR, nav_vol_tmp(:, idxn, ll,cc,jj,ii,kk).*exp(-1i*(read_shift.*kxLR+phase_shift.*kyLR+2/(2*LR_SIZE_Z)*kzLR(idxn))*2*pi), wwLR, LR_SIZE_XY, 2.5);
                    end
                    nav_imgb(:,:,:,cc) = fftshift(ifft(ifftshift(t_nav, 3), LR_SIZE_Z, 3),3);
                end
                
                
                dPhi(:) = dPhi(:) - dPhi0(:);
                dkx(:) = dkx(:) - dkx0(:);
                dky(:) = dky(:) - dky0(:);
                dkz(:) = dkz(:) - dkz0(:);
                
                % Take the mean of the initialization
                init_dkx = mean(col(dkx(:)));
                init_dky = mean(col(dky(:)));
                init_dkz = mean(col(dkz(:)));
                d_phi_tmp = sum(exp(1j*dPhi(:)));
                if abs(d_phi_tmp) < 0.1
                    d_phi_tmp =1;
                    fprintf(1,' *Warning: changing initialization phase to zero')
                end
                init_dPhi = angle(d_phi_tmp);
                init_parm = [2*pi*init_dky 2*pi*init_dkx 2*pi*init_dkz init_dPhi];
                
                
                f = @(s_vec)lsqcost_complex_mchan_Joe(s_vec, nav_imgb(:,:,:,:), nav_imgb0(:,:,:,:),mask);
                
                [x,resnorm,residual,exitflag,output,lambda,jacobian] =lsqnonlin(f, init_parm);
                s_vec = col(x);
                
                % Results go here
                corrections(1,ll,jj,ii,kk) = s_vec(2)/2/pi; %skx;
                corrections(2,ll,jj,ii,kk) = s_vec(1)/2/pi;% sky;
                corrections(3,ll,jj,ii,kk) = s_vec(3)/2/pi;% skz;
                corrections(4,ll,jj,ii,kk) = s_vec(4); %cPhi;
                
                clear s_vec;
                
            end
        end
    end
end

save recon_info/kspace_Corrections_meanshot.mat corrections

corrections(3,:,:,:,:) = corrections(3,:,:,:,:)*(-1);
save recon_info/kspace_Corrections_meanshot_flipz.mat corrections
        
        

       