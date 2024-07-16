Muzmat = load('SI/Mu.mat');
Muz = Muzmat.Mu;
Muxmat = load('LR/Mu.mat');
Mux = Muxmat.Mu;
Muymat = load('AP/Mu.mat');
Muy = Muymat.Mu;

XYdiff = Mux - Muy;
XZdiff = Mux - Muz;
YZdiff = Muy - Muz;

load('AP.mat', 'Zmotion')
load('AP.mat', 'Xmotion')
load('AP.mat', 'Ymotion')
YZmotion = Zmotion;
YXmotion = Xmotion;
YYmotion = Ymotion;
load('LR.mat', 'Zmotion')
load('LR.mat', 'Xmotion')
load('LR.mat', 'Ymotion')
XZmotion = Zmotion;
XXmotion = Xmotion;
XYmotion = Ymotion;
load('SI.mat', 'Zmotion')
load('SI.mat', 'Xmotion')
load('SI.mat', 'Ymotion')
ZZmotion = Zmotion;
ZXmotion = Xmotion;
ZYmotion = Ymotion;

