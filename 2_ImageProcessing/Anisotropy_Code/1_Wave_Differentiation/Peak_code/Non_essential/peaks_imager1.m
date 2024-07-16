clear all
dfile = dir('Multiexcit*.mat');
dfile2 = dir('BMR.mat');
dfile3 = dir('BuckyOut.mat');
load(dfile.name)
load(dfile2.name)
load(dfile3.name)

%peak1 = cat(3,peak1_X,peak1_Y,peak1_Z);
%peak2 = cat(3,peak2_X,peak2_Y,peak2_Z);
NUMROW = 80;
NUMCOL = 80;

peak1_R = A1x;
peak2_R = A2x;

peak_sum = peak1_R+peak2_R;
peak1_R_norm = peak1_R./peak_sum;
peak2_R_norm = peak2_R./peak_sum;

peak1_X = V1x(:,:,31,1);
peak1_Y = V1x(:,:,31,2);
peak1_Z = V1x(:,:,31,3);
peak2_X = V2x(:,:,31,1);
peak2_Y = V2x(:,:,31,2);
peak2_Z = V2x(:,:,31,3);

NNU(:,:,1,1) = peak1_X;
NNU(:,:,1,2) = peak1_Y;
NNU(:,:,1,3) = peak1_Z;
NNU(:,:,2,1) = peak2_X;
NNU(:,:,2,2) = peak2_Y;
NNU(:,:,2,3) = peak2_Z;
for ii = 1:NUMROW
for jj = 1:NUMCOL
    for kk = 1:2
NNUnorm(ii,jj,kk) = norm(squeeze(NNU(ii,jj,kk,:)));
    end
end
end

nscal = 2;
% SET UP FOR PLOTTING
% for color maps
% Propagation direction vectors
% FOR RGB COLOR MAPS OF PROPAGATION DIRECTION
imn(:,:,1) = nscal*(NNU(:,:,1,1)./NNUnorm(:,:,1));
imn(:,:,2) = nscal*(NNU(:,:,1,2)./NNUnorm(:,:,1));
imn(:,:,3) = nscal*(NNU(:,:,1,3)./NNUnorm(:,:,1));
figure;image(abs(imn)),axis equal,axis off

nscal = 2;
imn(:,:,1) = nscal*(NNU(:,:,2,1)./NNUnorm(:,:,2));
imn(:,:,2) = nscal*(NNU(:,:,2,2)./NNUnorm(:,:,2));
imn(:,:,3) = nscal*(NNU(:,:,2,3)./NNUnorm(:,:,2));
figure;image(abs(imn)),axis equal,axis off

