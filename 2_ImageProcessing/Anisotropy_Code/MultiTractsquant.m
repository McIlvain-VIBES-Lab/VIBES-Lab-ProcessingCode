tmp = load_nii('JHU-ICBM-labels-1mm.nii');
tmp = tmp.img;
img = zeros(size(tmp)); GenuCC = zeros(size(tmp)); BodyCC = zeros(size(tmp)); SpleniumCC = zeros(size(tmp));
CSTR = zeros(size(tmp)); CSTL = zeros(size(tmp)); ACRR = zeros(size(tmp)); ACRL = zeros(size(tmp));
SCRR = zeros(size(tmp)); SCRL = zeros(size(tmp)); PCRR = zeros(size(tmp)); PCRL = zeros(size(tmp));
SLFR = zeros(size(tmp)); SLFL = zeros(size(tmp));

for i = 1:48
    for l = 1:length(tmp(:,1,1))
        for j = 1:length(tmp(1,:,1))
            for k = 1:length(tmp(1,1,:))
                if tmp(l,j,k) == i
                    img(l,j,k,i) = 1;
                end
                if tmp(l,j,k) == 3
                    GenuCC(l,j,k) = 1;
                elseif tmp(l,j,k) == 4
                    BodyCC(l,j,k) = 1;
                elseif tmp(l,j,k) == 5
                    SpleniumCC(l,j,k) = 1;
                elseif tmp(l,j,k) == 7
                    CSTR(l,j,k) = 1;
                elseif tmp(l,j,k) == 8
                    CSTL(l,j,k) = 1;
                elseif tmp(l,j,k) == 23
                    ACRR(l,j,k) = 1;
                elseif tmp(l,j,k) == 24
                    ACRL(l,j,k) = 1;
                elseif tmp(l,j,k) == 25
                    SCRR(l,j,k) = 1;
                elseif tmp(l,j,k) == 26
                    SCRL(l,j,k) = 1;
                elseif tmp(l,j,k) == 27
                    PCRR(l,j,k) = 1;
                elseif tmp(l,j,k) == 28
                    PCRL(l,j,k) = 1;
                elseif tmp(l,j,k) == 41
                    SLFR(l,j,k) = 1;
                elseif tmp(l,j,k) == 42
                    SLFL(l,j,k) = 1;
                end
            end
        end
    end
end

%% Make .mat files into .nii

tGenuCC = make_nii(GenuCC); save_nii(tGenuCC,'tGenuCC.nii')
tBodyCC = make_nii(BodyCC); save_nii(tBodyCC,'tBodyCC.nii')
tSpleniumCC = make_nii(SpleniumCC); save_nii(tSpleniumCC,'tSpleniumCC.nii')
tCSTR = make_nii(CSTR); save_nii(tCSTR,'tCSTR.nii')
tCSTL = make_nii(CSTL); save_nii(tCSTL,'tCSTL.nii')
tACRR = make_nii(ACRR); save_nii(tACRR,'tACRR.nii')
tACRL = make_nii(ACRL); save_nii(tACRL,'tACRL.nii')
tSCRR = make_nii(SCRR); save_nii(tSCRR,'tSCRR.nii')
tSCRL = make_nii(SCRL); save_nii(tSCRL,'tSCRL.nii')
tPCRR = make_nii(PCRR); save_nii(tPCRR,'tPCRR.nii')
tPCRL = make_nii(PCRL); save_nii(tPCRL,'tPCRL.nii')
tSLFR = make_nii(SLFR); save_nii(tSLFR,'tSLFR.nii')
tSLFL = make_nii(SLFL); save_nii(tSLFL,'tSLFL.nii')

%% Flirting

!$FSLDIR/bin/flirt -in MNI152_T1_1mm_brain.nii -ref t2bet.nii -o ref.nii -omat output.txt
!$FSLDIR/bin/flirt -in tGenuCC.nii -ref ref.nii -o GenuCC.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tBodyCC.nii -ref ref.nii -o BodyCC.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tSpleniumCC.nii -ref ref.nii -o SpleniumCC.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tCSTR.nii -ref ref.nii -o CSTR.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tCSTL.nii -ref ref.nii -o CSTL.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tACRR.nii -ref ref.nii -o ACRR.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tACRL.nii -ref ref.nii -o ACRL.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tSCRR.nii -ref ref.nii -o SCRR.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tSCRL.nii -ref ref.nii -o SCRL.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tPCRR.nii -ref ref.nii -o PCRR.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tPCRL.nii -ref ref.nii -o PCRL.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tSLFR.nii -ref ref.nii -o SLFR.nii -init output.txt -applyxfm
!$FSLDIR/bin/flirt -in tSLFL.nii -ref ref.nii -o SLFL.nii -init output.txt -applyxfm
!gunzip *.nii.gz

%% Convert back from .nii to .mat

GenuCC = load_untouch_nii('GenuCC.nii'); GenuCC = GenuCC.img>0; BodyCC = load_untouch_nii('BodyCC.nii'); BodyCC = BodyCC.img>0;
SpleniumCC = load_untouch_nii('SpleniumCC.nii'); SpleniumCC = SpleniumCC.img>0;
CSTR = load_untouch_nii('CSTR.nii'); CSTR = CSTR.img>0; CSTL = load_untouch_nii('CSTL.nii'); CSTL = CSTL.img>0;
ACRR = load_untouch_nii('ACRR.nii'); ACRR = ACRR.img>0; ACRL = load_untouch_nii('ACRL.nii'); ACRL = ACRL.img>0;
SCRR = load_untouch_nii('SCRR.nii'); SCRR = SCRR.img>0; SCRL = load_untouch_nii('SCRL.nii'); SCRL = SCRL.img>0;
PCRR = load_untouch_nii('PCRR.nii'); PCRR = PCRR.img>0; PCRL = load_untouch_nii('PCRL.nii'); PCRL = PCRL.img>0;
SLFR = load_untouch_nii('SLFR.nii'); SLFR = SLFR.img>0; SLFL = load_untouch_nii('SLFL.nii'); SLFL = SLFL.img>0;

%% Move files
mkdir('flirtTracts')
mkdir('niiTracts')
movefile('tGenuCC.nii','flirtTracts/tGenuCC.nii'); movefile('GenuCC.nii','niiTracts/GenuCC.nii'); 
movefile('tBodyCC.nii','flirtTracts/tBodyCC.nii'); movefile('BodyCC.nii','niiTracts/BodyCC.nii'); 
movefile('tSpleniumCC.nii','flirtTracts/tSpleniumCC.nii'); movefile('SpleniumCC.nii','niiTracts/SpleniumCC.nii'); 
movefile('tCSTR.nii','flirtTracts/tCSTR.nii'); movefile('CSTR.nii','niiTracts/CSTR.nii'); 
movefile('tCSTL.nii','flirtTracts/tCSTL.nii'); movefile('CSTL.nii','niiTracts/CSTL.nii'); 
movefile('tACRR.nii','flirtTracts/tACRR.nii'); movefile('ACRR.nii','niiTracts/ACRR.nii'); 
movefile('tACRL.nii','flirtTracts/tACRL.nii'); movefile('ACRL.nii','niiTracts/ACRL.nii'); 
movefile('tSCRR.nii','flirtTracts/tSCRR.nii'); movefile('SCRR.nii','niiTracts/SCRR.nii'); 
movefile('tSCRL.nii','flirtTracts/tSCRL.nii'); movefile('SCRL.nii','niiTracts/SCRL.nii'); 
movefile('tPCRR.nii','flirtTracts/tPCRR.nii'); movefile('PCRR.nii','niiTracts/PCRR.nii'); 
movefile('tPCRL.nii','flirtTracts/tPCRL.nii'); movefile('PCRL.nii','niiTracts/PCRL.nii'); 
movefile('tSLFR.nii','flirtTracts/tSLFR.nii'); movefile('SLFR.nii','niiTracts/SLFR.nii'); 
movefile('tSLFL.nii','flirtTracts/tSLFL.nii'); movefile('SLFL.nii','niiTracts/SLFL.nii'); 
cd niiTracts
niis = dir('*.nii');

Tracts = cat(4,ACRL,ACRR,BodyCC,CSTL,CSTR,GenuCC,PCRL,PCRR,SCRL,SCRR,SLFL,SLFR,SpleniumCC);
cd ..
save Tracts.mat Tracts

cd ..; cd ..
load dtiall.mat
DTI = abs(dtiall(:,:,:,1))>0;

%% TRACTS
% 1 = Anterior Corona Radiata Left
% 2 = Anterior Corona Radiata Right
% 3 = Body Corpus Callosum
% 4 = Corticolspinal Tract Left
% 5 = Corticolspinal Tract Right
% 6 = Genu Corpus Callosum
% 7 = Posterior Corona Radiata Left
% 8 = Posterior Corona Radiata Right
% 9 = Superior Corona Radiata Left
% 10 = Superior Corona Radiata Right
% 11 = Superior Longitudinal Fasciculus Left
% 12 = Superior Longitudinal Fasciculus Right
% 13 = Splenium Corpus Callosum


cd Tracts/WMminor; load Tracts.mat
cd niiTracts; niis = dirs('*.nii');
cd ..; cd ..; cd ..
load dtiall.mat
DTI = abs(dtiall(:,:,:,1))>0;

%% 
load FA.mat
cd MRE; cd AP 
addpath('PostInv')
load Mu.mat; load slowfasttwo.mat; load BuckyOut.mat;
MutotlistAP = zeros(1,13); TotalAP = zeros(13,3); AmplistAP = zeros(1,13,2);
UperctotAP = zeros(1,13,4); ThlistAP = zeros(1,2,13); FAtractlistAP = zeros(1,13);


Theta1 = acos(dot(dtiall,flip(flip(permute(V1x,[2 1 3 4]),1),2),4));
Theta2 = acos(dot(dtiall,flip(flip(permute(V2x,[2 1 3 4]),1),2),4));

for ii = 1:length(niis)
    nii = niis(ii).name;
    tract = nii(1:end-4);
    TR = double(Tracts(:,:,:,ii));
    maskx = TR.*DTI;
    list = maskx(maskx==1);
    
    % Mu
    Mux = flip(flip(permute(Mu,[2 1 3]),1),2);
    nn = Mux.*maskx; Mulist = nn(nn>0);
    MutotlistAP(1:length(Mulist),ii) = Mulist;
    Mumean = mean(Mulist);
    
    %Theta
    Theta1mask = Theta1.*maskx; Theta1list = Theta1(Theta1mask>0);
    Theta2mask = Theta2.*maskx; Theta2list = Theta2(Theta2mask>0);
    
    ThlistAP(1:length(Theta1list),1,ii) = Theta1list;
    ThlistAP(1:length(Theta1list),2,ii) = Theta2list;
    
    %Amplitude weightings
    A1p = A1x./(A1x+A2x); A2p = A2x./(A1x+A2x);
    A1mask = A1p.*maskx; A1list = A1mask(A1mask>0);
    A2mask = A2p.*maskx; A2list = A2mask(A2mask>0);
    
    AmplistAP(1:length(A1list),1,ii) = A1list;
    AmplistAP(1:length(A1list),2,ii) = A2list;
    
    %FA
    FAmask = FA.*maskx; FAlist = FA(FAmask>0);
    FAtractlistAP(1:length(FAlist),ii) = FAlist;
    
    %Slow/Fast
    UsV1 = flip(flip(permute(U_s_V1,[2 1 3]),1),2); UfV1 = flip(flip(permute(U_f_V1,[2 1 3]),1),2);
    UsV2 = flip(flip(permute(U_s_V2,[2 1 3]),1),2); UfV2 = flip(flip(permute(U_f_V2,[2 1 3]),1),2);
    
    UsV1list = abs(UsV1(maskx>0)); UsV1list = UsV1list(~isnan(UsV1list));
    UfV1list = abs(UfV1(maskx>0)); UfV1list = UfV1list(~isnan(UfV1list));
    UsV2list = abs(UsV2(maskx>0)); UsV2list = UsV2list(~isnan(UsV2list));
    UfV2list = abs(UfV2(maskx>0)); UfV2list = UfV2list(~isnan(UfV2list));
    
    UspercV1 = UsV1list./(UfV1list+UsV1list)*100; UfpercV1 = UfV1list./(UfV1list+UsV1list)*100;
    UspercV2 = UsV2list./(UfV2list+UsV2list)*100; UfpercV2 = UfV2list./(UfV2list+UsV2list)*100;
    
    UperctotAP(1:length(UspercV1),ii,3) = UspercV2; UperctotAP(1:length(UspercV1),ii,4) = UfpercV2;
    UperctotAP(1:length(UspercV1),ii,1) = UspercV1; UperctotAP(1:length(UspercV1),ii,2) = UfpercV1;
    
    Us = flip(flip(permute(U_s>U_f,[2 1 3]),1),2); Uf = flip(flip(permute(U_f>U_s,[2 1 3]),1),2);
    Usmask = maskx.*Us; Ufmask = maskx.*Uf; masktot = Usmask+Ufmask;
    Slow = Us(Usmask == 1); Slowperc = length(Slow)/length(maskx(maskx==1))*100;
    Fast = Uf(Ufmask == 1); Fastperc = length(Fast)/length(maskx(maskx==1))*100;
    
    TotalAP(ii,1) = Mumean; TotalAP(ii,2) = Slowperc; TotalAP(ii,3) = Fastperc;
    
    disp(sprintf('Mulist: %i; UsV1list: %i; Theta1list: %i; ii: %i',length(Mulist),length(UsV1list),length(Theta1list),ii))
    
end

save('Quant_MUFS_AP.mat','TotalAP','UperctotAP','MutotlistAP','ThlistAP','FAtractlistAP')

cd ..; cd LR
addpath('PostInv')
load Mu.mat; load slowfasttwo.mat; load BuckyOut.mat
MutotlistLR = zeros(1,13); TotalLR = zeros(13,3); AmplistLR = zeros(1,13,2);
UperctotLR = zeros(1,13,4); ThlistLR = zeros(1,2,13); FAtractlistLR = zeros(1,13);

Theta1 = acos(dot(dtiall,flip(flip(permute(V1x,[2 1 3 4]),1),2),4));
Theta2 = acos(dot(dtiall,flip(flip(permute(V2x,[2 1 3 4]),1),2),4));

for ii = 1:length(niis)
    nii = niis(ii).name;
    tract = nii(1:end-4);
    TR = double(Tracts(:,:,:,ii));
    maskx = TR.*DTI;
    list = maskx(maskx==1);
    
    % Mu
    Mux = flip(flip(permute(Mu,[2 1 3]),1),2);
    nn = Mux.*maskx; Mulist = nn(nn>0);
    MutotlistLR(1:length(Mulist),ii) = Mulist;
    Mumean = mean(Mulist);
    
    %Theta
    Theta1mask = Theta1.*maskx; Theta1list = Theta1(Theta1mask>0);
    Theta2mask = Theta2.*maskx; Theta2list = Theta2(Theta2mask>0);
    
    ThlistLR(1:length(Theta1list),1,ii) = Theta1list;
    ThlistLR(1:length(Theta1list),2,ii) = Theta2list;
    
    %Amplitude weightings
    A1p = A1x./(A1x+A2x); A2p = A2x./(A1x+A2x);
    A1mask = A1p.*maskx; A1list = A1mask(A1mask>0);
    A2mask = A2p.*maskx; A2list = A2mask(A2mask>0);
    
    AmplistLR(1:length(A1list),1,ii) = A1list;
    AmplistLR(1:length(A1list),2,ii) = A2list;
    
    %FA
    FAmask = FA.*maskx; FAlist = FA(FAmask>0);
    FAtractlistLR(1:length(FAlist),ii) = FAlist;
    
    %Slow/Fast
    UsV1 = flip(flip(permute(U_s_V1,[2 1 3]),1),2); UfV1 = flip(flip(permute(U_f_V1,[2 1 3]),1),2);
    UsV2 = flip(flip(permute(U_s_V2,[2 1 3]),1),2); UfV2 = flip(flip(permute(U_f_V2,[2 1 3]),1),2);
    
    UsV1list = abs(UsV1(maskx>0)); UsV1list = UsV1list(~isnan(UsV1list));
    UfV1list = abs(UfV1(maskx>0)); UfV1list = UfV1list(~isnan(UfV1list));
    UsV2list = abs(UsV2(maskx>0)); UsV2list = UsV2list(~isnan(UsV2list));
    UfV2list = abs(UfV2(maskx>0)); UfV2list = UfV2list(~isnan(UfV2list));
    
    UspercV1 = UsV1list./(UfV1list+UsV1list)*100; UfpercV1 = UfV1list./(UfV1list+UsV1list)*100;
    UspercV2 = UsV2list./(UfV2list+UsV2list)*100; UfpercV2 = UfV2list./(UfV2list+UsV2list)*100;
    
    UperctotLR(1:length(UspercV1),ii,3) = UspercV2; UperctotLR(1:length(UspercV1),ii,4) = UfpercV2;
    UperctotLR(1:length(UspercV1),ii,1) = UspercV1; UperctotLR(1:length(UspercV1),ii,2) = UfpercV1;
    
    Us = flip(flip(permute(U_s>U_f,[2 1 3]),1),2); Uf = flip(flip(permute(U_f>U_s,[2 1 3]),1),2);
    Usmask = maskx.*Us; Ufmask = maskx.*Uf; masktot = Usmask+Ufmask;
    Slow = Us(Usmask == 1); Slowperc = length(Slow)/length(maskx(maskx==1))*100;
    Fast = Uf(Ufmask == 1); Fastperc = length(Fast)/length(maskx(maskx==1))*100;
    
    TotalLR(ii,1) = Mumean; TotalLR(ii,2) = Slowperc; TotalLR(ii,3) = Fastperc;
    
    disp(sprintf('Mulist: %i; UsV1list: %i; Theta1list: %i; ii: %i',length(Mulist),length(UsV1list),length(Theta1list),ii))
    
end

save('Quant_MUFS_LR.mat','TotalLR','UperctotLR','MutotlistLR','ThlistLR','FAtractlistLR')

cd ..; cd SI
addpath('PostInv')
load Mu.mat; load slowfasttwo.mat; load BuckyOut.mat
MutotlistSI = zeros(1,13); TotalSI = zeros(13,3); AmplistSI = zeros(1,13,2);
UperctotSI = zeros(1,13,4); ThlistSI = zeros(1,2,13); FAtractlistSI = zeros(1,13);

Theta1 = acos(dot(dtiall,flip(flip(permute(V1x,[2 1 3 4]),1),2),4));
Theta2 = acos(dot(dtiall,flip(flip(permute(V2x,[2 1 3 4]),1),2),4));

for ii = 1:length(niis)
    nii = niis(ii).name;
    tract = nii(1:end-4);
    TR = double(Tracts(:,:,:,ii));
    maskx = TR.*DTI;
    list = maskx(maskx==1);
    
    % Mu
    Mux = flip(flip(permute(Mu,[2 1 3]),1),2);
    nn = Mux.*maskx; Mulist = nn(nn>0);
    MutotlistSI(1:length(Mulist),ii) = Mulist;
    Mumean = mean(Mulist);
    
    %Theta
    Theta1mask = Theta1.*maskx; Theta1list = Theta1(Theta1mask>0);
    Theta2mask = Theta2.*maskx; Theta2list = Theta2(Theta2mask>0);
    
    ThlistSI(1:length(Theta1list),1,ii) = Theta1list;
    ThlistSI(1:length(Theta1list),2,ii) = Theta2list;
  
    %Amplitude weightings
    A1p = A1x./(A1x+A2x); A2p = A2x./(A1x+A2x);
    A1mask = A1p.*maskx; A1list = A1mask(A1mask>0);
    A2mask = A2p.*maskx; A2list = A2mask(A2mask>0);
    
    AmplistSI(1:length(A1list),1,ii) = A1list;
    AmplistSI(1:length(A1list),2,ii) = A2list;
    
    %FA
    FAmask = FA.*maskx; FAlist = FA(FAmask>0);
    FAtractlistSI(1:length(FAlist),ii) = FAlist;
    
    %Slow/Fast
    UsV1 = flip(flip(permute(U_s_V1,[2 1 3]),1),2); UfV1 = flip(flip(permute(U_f_V1,[2 1 3]),1),2);
    UsV2 = flip(flip(permute(U_s_V2,[2 1 3]),1),2); UfV2 = flip(flip(permute(U_f_V2,[2 1 3]),1),2);
    
    UsV1list = abs(UsV1(maskx>0)); UsV1list = UsV1list(~isnan(UsV1list));
    UfV1list = abs(UfV1(maskx>0)); UfV1list = UfV1list(~isnan(UfV1list));
    UsV2list = abs(UsV2(maskx>0)); UsV2list = UsV2list(~isnan(UsV2list));
    UfV2list = abs(UfV2(maskx>0)); UfV2list = UfV2list(~isnan(UfV2list));
    
    UspercV1 = UsV1list./(UfV1list+UsV1list)*100; UfpercV1 = UfV1list./(UfV1list+UsV1list)*100;
    UspercV2 = UsV2list./(UfV2list+UsV2list)*100; UfpercV2 = UfV2list./(UfV2list+UsV2list)*100;
    
    UperctotSI(1:length(UspercV1),ii,3) = UspercV2; UperctotSI(1:length(UspercV1),ii,4) = UfpercV2;
    UperctotSI(1:length(UspercV1),ii,1) = UspercV1; UperctotSI(1:length(UspercV1),ii,2) = UfpercV1;
    
    Us = flip(flip(permute(U_s>U_f,[2 1 3]),1),2); Uf = flip(flip(permute(U_f>U_s,[2 1 3]),1),2);
    Usmask = maskx.*Us; Ufmask = maskx.*Uf; masktot = Usmask+Ufmask;
    Slow = Us(Usmask == 1); Slowperc = length(Slow)/length(maskx(maskx==1))*100;
    Fast = Uf(Ufmask == 1); Fastperc = length(Fast)/length(maskx(maskx==1))*100;
    
    TotalSI(ii,1) = Mumean; TotalSI(ii,2) = Slowperc; TotalSI(ii,3) = Fastperc;
    
    disp(sprintf('Mulist: %i; UsV1list: %i; Theta1list: %i; ii: %i',length(Mulist),length(UsV1list),length(Theta1list),ii))
    
end

save('Quant_MUFS_SI.mat','TotalSI','UperctotSI','MutotlistSI','ThlistSI','FAtractlistSI')

