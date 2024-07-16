% three excitations, slow/fast discrimination on V1
% discussed by DRS and CLJ, 1/9/2018, we think this code is good
for ti = 1:13;
SF_thresh = 70;
FA_thresh = 0.7;

FAv = FAtractlistAP(:,ti)>FA_thresh;

Sv1_AP = UperctotAP(:,ti,1)>=SF_thresh;
Fv1_AP = UperctotAP(:,ti,2)>=SF_thresh;
MuS_AP = MutotlistAP(logical(Sv1_AP.*FAv),ti);
MuF_AP = MutotlistAP(logical(Fv1_AP.*FAv),ti);

Sv1_LR = UperctotLR(:,ti,1)>=SF_thresh;
Fv1_LR = UperctotLR(:,ti,2)>=SF_thresh;
MuS_LR = MutotlistLR(logical(Sv1_LR.*FAv),ti);
MuF_LR = MutotlistLR(logical(Fv1_LR.*FAv),ti);

Sv1_SI = UperctotSI(:,ti,1)>=SF_thresh;
Fv1_SI = UperctotSI(:,ti,2)>=SF_thresh;
MuS_SI = MutotlistSI(logical(Sv1_SI.*FAv),ti);
MuF_SI = MutotlistSI(logical(Fv1_SI.*FAv),ti);

b = [MuS_LR;MuF_LR;MuS_AP;MuF_AP;MuS_SI;MuF_SI];

AS_LR(1:length(MuS_LR),1) = 1;
AS_LR(1:length(MuS_LR),2) = cos(ThlistLR(logical(Sv1_LR.*FAv),1,ti)).^2;
AS_LR(1:length(MuS_LR),3) = 0;

AF_LR(1:length(MuF_LR),1) = 1;
AF_LR(1:length(MuF_LR),2) = cos(2*ThlistLR(logical(Fv1_LR.*FAv),1,ti)).^2;
AF_LR(1:length(MuF_LR),3) = sin(2*ThlistLR(logical(Fv1_LR.*FAv),1,ti)).^2;

AS_AP(1:length(MuS_AP),1) = 1;
AS_AP(1:length(MuS_AP),2) = cos(ThlistAP(logical(Sv1_AP.*FAv),1,ti)).^2;
AS_AP(1:length(MuS_AP),3) = 0;

AF_AP(1:length(MuF_AP),1) = 1;
AF_AP(1:length(MuF_AP),2) = cos(2*ThlistAP(logical(Fv1_AP.*FAv),1,ti)).^2;
AF_AP(1:length(MuF_AP),3) = sin(2*ThlistAP(logical(Fv1_AP.*FAv),1,ti)).^2;

AS_SI(1:length(MuS_SI),1) = 1;
AS_SI(1:length(MuS_SI),2) = cos(ThlistSI(logical(Sv1_SI.*FAv),1,ti)).^2;
AS_SI(1:length(MuS_SI),3) = 0;

AF_SI(1:length(MuF_SI),1) = 1;
AF_SI(1:length(MuF_SI),2) = cos(2*ThlistSI(logical(Fv1_SI.*FAv),1,ti)).^2;
AF_SI(1:length(MuF_SI),3) = sin(2*ThlistSI(logical(Fv1_SI.*FAv),1,ti)).^2;

A = cat(1,AS_LR,AF_LR,AS_AP,AF_AP,AS_SI,AF_SI);
x = A\b;
b_new = A*x;
x_soln(:,ti) = x;
r2(ti,1) = rsquare(b,b_new);
clear b x A b_new Sv1_LR Fv1_LR Sv1_AP Fv1_AP AS_LR AF_LR AS_AP AF_AP MuS_LR MuF_LR MuS_AP MuF_AP Sv1_SI Fv1_SI AS_SI AF_SI MuS_SI MuF_SI
end