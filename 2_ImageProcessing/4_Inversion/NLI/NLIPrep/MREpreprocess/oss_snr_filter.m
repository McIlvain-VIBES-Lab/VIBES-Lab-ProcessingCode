function OSS_SNR = oss_snr_filter(pimg,voxres,freq,mask)

[ny,nx,nz,ni,nj] = size(pimg);

phs_img = zeros(size(pimg));
for jj = 1:nj
    for ii = 1:ni
        tmp1 = pimg(:,:,:,ii,jj);
%         tmp2 = col(tmp1(mask>0));
        tmp2=tmp1(mask>0);
        tmp3 = fft(tmp2);
        tmp3(1) = 0;
        tmp4 = ifft(tmp3);
        phs_img(:,:,:,ii,jj) = embed(tmp4,logical(mask));
    end
end

fft_img = fft(phs_img,[],5)/(nj/2);
wave_img = fft_img(:,:,:,:,2);

t= (0:2*pi/nj:(2*nj-1)/nj*pi)';
tt = repmat(reshape(t,[1 1 1 1 nj]),[ny nx nz ni 1]);
err_img = std(real(repmat(wave_img,[1 1 1 1 nj]).*exp(1i*tt)) - phs_img,0,5);

convdata(voxres,freq,wave_img,err_img,mask,mask);
cd dcfiles
if nj == 8
    [OSS_SNR,~,~,~,~,~]=Strain_SNR;
    OSS_SNR = OSS_SNR*1.35;
elseif nj == 4
    [OSS_SNR,~,~,~,~,~]=Strain_SNR;
    OSS_SNR = OSS_SNR*1.35/(sqrt(2.5));
end
cd ..