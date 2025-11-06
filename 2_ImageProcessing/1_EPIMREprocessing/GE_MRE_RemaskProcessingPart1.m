function [phsimg] = GE_MRE_RemaskProcessingPart1(dx,dy,dz,freq)
%% Section 1 of the MRE Processing Code
% GE_MRE_Processing Part 1
% March 26th 2025
% Grace McIlvain

% This code pulls all dicoms and makes t2stack


dirlist2 = dir('IM*');
nn = 0;
for bb =1:4
for aa = 1:(length(dirlist2)/4) %% number of slices
%    for dd = 1:45
        nn=nn+1;
        images1(:,:,aa,bb) = dicomread(dirlist2(nn).name);
    end
end

tt = 1;
rr = 1;
for jj=1:(length(dirlist2)/48)
for ii=1:4

e1_r(:,:,rr,ii) = images1(:,:,tt,ii);
e1_i(:,:,rr,ii) = images1(:,:,tt+1,ii);

e2_r(:,:,rr,ii) = images1(:,:,tt+2,ii);
e2_i(:,:,rr,ii) = images1(:,:,tt+3,ii);

e3_r(:,:,rr,ii) = images1(:,:,tt+4,ii);
e3_i(:,:,rr,ii) = images1(:,:,tt+5,ii);

e4_r(:,:,rr,ii) = images1(:,:,tt+6,ii);
e4_i(:,:,rr,ii) = images1(:,:,tt+7,ii);

e5_r(:,:,rr,ii) = images1(:,:,tt+8,ii);
e5_i(:,:,rr,ii) = images1(:,:,tt+9,ii);

e6_r(:,:,rr,ii) = images1(:,:,tt+10,ii);
e6_i(:,:,rr,ii) = images1(:,:,tt+11,ii);

end
tt = tt+12;
rr =rr+1;
end


e1_t = double(e1_r)+i*double(e1_i);
e2_t = double(e2_r)+i*double(e2_i);
e3_t = double(e3_r)+i*double(e3_i);
e4_t = double(e4_r)+i*double(e4_i);
e5_t = double(e5_r)+i*double(e5_i);
e6_t = double(e6_r)+i*double(e6_i);

e1 = flip(flip(permute(e1_t,[2 1 3 4]),2),1);
e2 = flip(flip(permute(e2_t,[2 1 3 4]),2),1);
e3 = flip(flip(permute(e3_t,[2 1 3 4]),2),1);
e4 = flip(flip(permute(e4_t,[2 1 3 4]),2),1);
e5 = flip(flip(permute(e5_t,[2 1 3 4]),2),1);
e6 = flip(flip(permute(e6_t,[2 1 3 4]),2),1);

cd ..

phsimg(:,:,:,3,:)=angle(e1./e2);
magimg(:,:,:,3,:)=((abs(e1)+abs(e2))/max(col(abs(e1)+abs(e2))))*30000;
phsimg(:,:,:,2,:)=angle(e3./e4);
magimg(:,:,:,2,:)=((abs(e3)+abs(e4))/max(col(abs(e3)+abs(e4))))*30000;
phsimg(:,:,:,1,:)=angle(e5./e6);
magimg(:,:,:,1,:)=((abs(e5)+abs(e6))/max(col(abs(e5)+abs(e6))))*30000;

% GM build in distortion correction
%GE_MRE_DistortionCorrection
%end

end