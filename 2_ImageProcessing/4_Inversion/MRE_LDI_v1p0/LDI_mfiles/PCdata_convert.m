tt = 0:(2*pi/8):(2*pi*(7/8))
Wdata = cat(5,Xmotion,Ymotion,Zmotion);
for ii = 1:length(tt)
PCdata(:,:,:,ii,:) = real(Wdata.*exp(-1i*tt(ii)));
end