function WMM_minor
%% Create minor Tract masks 

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
