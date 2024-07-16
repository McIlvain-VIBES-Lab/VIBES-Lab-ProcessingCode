% three excitations, slow/fast discrimination on V1
% discussed by DRS and CLJ, 1/9/2018, we think this code is good
   
SF_thresh = 60;
FA_thresh = 0.7;

FAv = permute(FAtractlistAP(1,:)>FA_thresh,[2 1]);

Sv1_AP = UperctotAP(:,1)>=SF_thresh;
Fv1_AP = UperctotAP(:,2)>=SF_thresh;
MuS_AP = permute(MutotlistAP(1,logical(Sv1_AP.*FAv)),[2 1]);
MuF_AP = permute(MutotlistAP(1,logical(Fv1_AP.*FAv)),[2 1]);

Sv1_LR = UperctotLR(:,1)>=SF_thresh;
Fv1_LR = UperctotLR(:,2)>=SF_thresh;
MuS_LR = permute(MutotlistLR(logical(Sv1_LR.*FAv)),[2 1]);
MuF_LR = permute(MutotlistLR(logical(Fv1_LR.*FAv)),[2 1]);

b = [MuS_LR;MuF_LR;MuS_AP;MuF_AP];

AS_LR(1:length(MuS_LR),1) = 1;
AS_LR(1:length(MuS_LR),2) = cos(ThlistLR(logical(Sv1_LR.*FAv),1)).^2;
AS_LR(1:length(MuS_LR),3) = 0;

AF_LR(1:length(MuF_LR),1) = 1;
AF_LR(1:length(MuF_LR),2) = cos(2*ThlistLR(logical(Fv1_LR.*FAv),1)).^2;
AF_LR(1:length(MuF_LR),3) = sin(2*ThlistLR(logical(Fv1_LR.*FAv),1)).^2;

AS_AP(1:length(MuS_AP),1) = 1;
AS_AP(1:length(MuS_AP),2) = cos(ThlistAP(logical(Sv1_AP.*FAv),1)).^2;
AS_AP(1:length(MuS_AP),3) = 0;

AF_AP(1:length(MuF_AP),1) = 1;
AF_AP(1:length(MuF_AP),2) = cos(2*ThlistAP(logical(Fv1_AP.*FAv),1)).^2;
AF_AP(1:length(MuF_AP),3) = sin(2*ThlistAP(logical(Fv1_AP.*FAv),1)).^2;

A = cat(1,AS_LR,AF_LR,AS_AP,AF_AP);
x = A\b;
b_new = A*x;
x_soln(:) = x;
r2(1) = rsquare(b,b_new);
clear b x A b_new Sv1_LR Fv1_LR Sv1_AP Fv1_AP AS_LR AF_LR AS_AP AF_AP MuS_LR MuF_LR MuS_AP MuF_AP Sv1_SI Fv1_SI AS_SI AF_SI MuS_SI MuF_SI