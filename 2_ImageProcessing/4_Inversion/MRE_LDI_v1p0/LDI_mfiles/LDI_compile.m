%% Load Data
load BuckyOut.mat
tmp = dir('Multi*');
tmp1 = tmp.name;
load(sprintf('%s',tmp1))

%% Mask data point collection
cd ..; cd ..; load dtiall.mat
dtimask = dtiall(:,:,:,1)>0;
cd Tracts/WMmain
CC = load_untouch_nii('tractsCC.nii');
CCmask = CC.img>0;
sesz = 2;
SE = strel('ball',sesz,sesz,6);
CCerode = imerode(double(CCmask),SE);
CCerode2 = (CCerode+sesz)>0.5;
se2 = strel('diamond',2);
CCerode3 = imerode(double(CCerode2),se2);
msk = dtimask.*CCerode3;
cd ..; cd ..;

cd MRE/AP
nn = 0;
clear lst
for ii = 1:length(msk(:,1,1))
    for jj = 1:length(msk(1,:,1))
        for kk = 1:length(msk(1,1,:))
            if msk(ii,jj,kk) == 1
                nn = nn + 1;
                lst(nn,:) = [ii jj kk];
            end
        end
    end
end

%% Stiffness map building for Slow 1 + 2 & Fast 1 + 2
Mu_s1 = zeros(size(msk)); Mu_f1 = zeros(size(msk));
Mu_s2 = zeros(size(msk)); Mu_f2 = zeros(size(msk));
Gps1x = zeros(size(msk)); Gpf1x = zeros(size(msk));
Gps2x = zeros(size(msk)); Gpf2x = zeros(size(msk));
Gdps1x = zeros(size(msk)); Gdpf1x = zeros(size(msk));
Gdps2x = zeros(size(msk)); Gdpf2x = zeros(size(msk));

for ii = 1:nn
    ii
    x = lst(ii,1); y = lst(ii,2); z = lst(ii,3);
    [Mus1,Muf1,Mus2,Muf2,Gp,Gdp] = skel_filter(x,y,z,V1index,V2index,mask,t2stack,dtiall);
    Mu_s1(x,y,z) = Mus1; 
    Mu_f1(x,y,z) = Muf1; 
    Mu_s2(x,y,z) = Mus2; 
    Mu_f2(x,y,z) = Muf2; 
    Gps1x(x,y,z) = Gp(1);
    Gpf1x(x,y,z) = Gp(2);
    Gps2x(x,y,z) = Gp(3);
    Gpf2x(x,y,z) = Gp(4);
    Gdps1x(x,y,z) = Gdp(1);
    Gdpf1x(x,y,z) = Gdp(2);
    Gdps2x(x,y,z) = Gdp(3);
    Gdpf2x(x,y,z) = Gdp(4);
end

save Mu_LDI_diff.mat Mu_s1 Mu_f1 Mu_s2 Mu_f2 Gps1x Gpf1x Gps2x Gpf2x Gdps1x Gdpf1x Gdps2x Gdpf2x
% 
% cd ..; cd LR
% 
% load BuckyOut.mat
% 
% nn = 0;
% clear lst
% for ii = 1:length(msk(:,1,1))
%     for jj = 1:length(msk(1,:,1))
%         for kk = 1:length(msk(1,1,:))
%             if msk(ii,jj,kk) == 1
%                 nn = nn + 1;
%                 lst(nn,:) = [ii jj kk];
%             end
%         end
%     end
% end
% 
% %% Stiffness map building for Slow 1 + 2 & Fast 1 + 2
% Mu_s1 = zeros(size(msk)); Mu_f1 = zeros(size(msk));
% Mu_s2 = zeros(size(msk)); Mu_f2 = zeros(size(msk));
% Gps1x = zeros(size(msk)); Gpf1x = zeros(size(msk));
% Gps2x = zeros(size(msk)); Gpf2x = zeros(size(msk));
% Gdps1x = zeros(size(msk)); Gdpf1x = zeros(size(msk));
% Gdps2x = zeros(size(msk)); Gdpf2x = zeros(size(msk));
% 
% for ii = 1:nn
%     ii
%     x = lst(ii,1); y = lst(ii,2); z = lst(ii,3);
%     [Mus1,Muf1,Mus2,Muf2,Gp,Gdp] = skel_filter1(x,y,z,V1index,V2index,mask,t2stack,dtiall);
%     Mu_s1(x,y,z) = Mus1; 
%     Mu_f1(x,y,z) = Muf1; 
%     Mu_s2(x,y,z) = Mus2; 
%     Mu_f2(x,y,z) = Muf2;
%     Gps1x(x,y,z) = Gp(1);
%     Gpf1x(x,y,z) = Gp(2);
%     Gps2x(x,y,z) = Gp(3);
%     Gpf2x(x,y,z) = Gp(4);
%     Gdps1x(x,y,z) = Gdp(1);
%     Gdpf1x(x,y,z) = Gdp(2);
%     Gdps2x(x,y,z) = Gdp(3);
%     Gdpf2x(x,y,z) = Gdp(4);
% end
% 
% save Mu_LDI_diff414.mat Mu_s1 Mu_f1 Mu_s2 Mu_f2 Gps1x Gpf1x Gps2x Gpf2x Gdps1x Gdpf1x Gdps2x Gdpf2x
