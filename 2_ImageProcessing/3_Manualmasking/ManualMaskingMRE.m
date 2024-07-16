%% Manual MRE image masking

load('t2mask_bet.mat')
load('t2stack.mat')
tmp = (t2stack.*mask)./(1000);
%tmp = (t2stack.*mask)./(2.8e-6); %Might have to change the denominator

maskx = zeros(size(mask));


% look at the whole brain
figure;im(tmp);caxis([0 1])

% look at the original mask
figure;im(mask)

% look at the mask you have created (the negative mask)
figure;im(maskx)

figure;im(tmp.*mask.*abs(1-maskx));caxis([0 1])


for ss = 23:-1:9 %Slice range
    ss
    maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
end
save maskx.mat maskx 




