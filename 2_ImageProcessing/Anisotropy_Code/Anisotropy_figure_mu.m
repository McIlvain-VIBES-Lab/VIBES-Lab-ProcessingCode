cd /Volumes/CLJ-001/Group/drsmitty/ScanData/Multi/170222/Multiexcite_01
[x_soln, crit_out1, x_soln1] = combo_2waves_LDI_crit(2,1,2);
mu1 = x_soln1(:,1,1);
phi1 = x_soln1(:,2,1);
zeta1 = x_soln1(:,3,1);
crit1 = crit_out1(:,:,1);
tmp1 = (mu1>1000).*(phi1>0).*(zeta1>0).*(phi1<(2*mu1)).*(zeta1<(2*mu1)).*(mu1<6000);
phi1 = phi1(logical(tmp1))./mu1(logical(tmp1));
zeta1 = zeta1(logical(tmp1))./mu1(logical(tmp1));
mu1 = mu1(logical(tmp1));
crit1 = crit1(logical(tmp1));

cd /Volumes/CLJ-001/Group/drsmitty/ScanData/Multi/170222/Multiexcite_02
[x_soln, crit_out2, x_soln2] = combo_2waves_LDI_crit(2,1,2);
mu2 = x_soln2(:,1,1);
phi2 = x_soln2(:,2,1);
zeta2 = x_soln2(:,3,1);
crit2 = crit_out2(:,:,1);
tmp2 = (mu2>1000).*(phi2>0).*(zeta2>0).*(phi2<(2*mu2)).*(zeta2<(2*mu2)).*(mu2<6000);
phi2 = phi2(logical(tmp2))./mu2(logical(tmp2));
zeta2 = zeta2(logical(tmp2))./mu2(logical(tmp2));
mu2 = mu2(logical(tmp2));
crit2 = crit2(logical(tmp2));


cd /Volumes/CLJ-001/Group/drsmitty/ScanData/Multi/170222/Multiexcite_03
[x_soln, crit_out3, x_soln3] = combo_2waves_LDI_crit(2,1,2);
mu3 = x_soln3(:,1,1);
phi3 = x_soln3(:,2,1);
zeta3 = x_soln3(:,3,1);
crit3 = crit_out3(:,:,1);
tmp3 = (mu3>1000).*(phi3>0).*(zeta3>0).*(phi3<(2*mu3)).*(zeta3<(2*mu3)).*(mu3<6000);
phi3 = phi3(logical(tmp3))./mu3(logical(tmp3));
zeta3 = zeta3(logical(tmp3))./mu3(logical(tmp3));
mu3 = mu3(logical(tmp3));
crit3 = crit3(logical(tmp3));

cd /Volumes/CLJ-001/Group/drsmitty/ScanData/Multi/170816/Clj_Multiexcite_04_02
[x_soln, crit_out4, x_soln4] = combo_2waves_LDI_crit(2,1,2);
mu4 = x_soln4(:,1,1);
phi4 = x_soln4(:,2,1);
zeta4 = x_soln4(:,3,1);
crit4 = crit_out4(:,:,1);
tmp4 = (mu4>1000).*(phi4>0).*(zeta4>0).*(phi4<(2*mu4)).*(zeta4<(2*mu4)).*(mu4<6000);
phi4 = phi4(logical(tmp4))./mu4(logical(tmp4));
zeta4 = zeta4(logical(tmp4))./mu4(logical(tmp4));
mu4 = mu4(logical(tmp4));
crit4 = crit4(logical(tmp4));

mu = cat(1,mu1,mu2,mu3,mu4);
phi = cat(1,phi1,phi2,phi3,phi4);
zeta = cat(1,zeta1,zeta2,zeta3,zeta4);
tmp = cat(1,tmp1,tmp2,tmp3,tmp4);
crit = cat(1,crit1,crit2,crit3,crit4);

% x_solna = cat(3,x_soln1,x_soln2,x_soln3,x_soln4);
% tmp = (mu1>0).*(phi1>0).*(zeta1>0);

median(mu)
median(zeta)
median(phi)
length(mu)
sum(tmp)

zetasort = sort(zeta)
phisort = sort(phi)
musort = sort(mu)


% critall = cat(1,crit_out1(:,:,1),crit_out2(:,:,1),crit_out3(:,:,1),crit_out4(:,:,1));

% muplus = mean(mu)+std(mu); phiplus = mean(phi)+std(phi); zetaplus = mean(zeta)+std(zeta);
% 
% munew = mu(mu<muplus); phinew = phi(phi<phiplus); zetanew = zeta(zeta<zetaplus);
% mumean = mean(munew); phimean = mean(phinew); zetamean = mean(zetanew);
% 
% mupn = mean(munew)+std(munew);
% phipn = mean(phinew)+std(munew);
% zetapn = mean(zetanew)+std(zetanew);
% 
% munn = munew(munew<mupn); phinn = phinew(phinew<phipn); zetann = zetanew(zetanew<zetapn);
% 
% mean(munn)
% mean(phinn)
% mean(zetann)

figure; boxplot([mu;mu],[ones(length(zeta1),1);2*ones(length(zeta2),1);3*ones(length(zeta3),1);4*ones(length(zeta4),1);5*ones(length(zeta),1)])
figure; boxplot([phi;phi],[ones(length(zeta1),1);2*ones(length(zeta2),1);3*ones(length(zeta3),1);4*ones(length(zeta4),1);5*ones(length(zeta),1)])
figure; boxplot([zeta;zeta],[ones(length(zeta1),1);2*ones(length(zeta2),1);3*ones(length(zeta3),1);4*ones(length(zeta4),1);5*ones(length(zeta),1)])

figure; plot(mu,phi)
figure; plot(mu,zeta)
figure; plot(phi,zeta)

figure;plot(phi1,zeta1)
hold on
plot(phi2,zeta2)
plot(phi3,zeta3)
plot(phi4,zeta4)


figure;notBoxPlot([mu;mu],[ones(length(zeta1),1);2*ones(length(zeta2),1);3*ones(length(zeta3),1);4*ones(length(zeta4),1);5*ones(length(zeta),1)])
figure;notBoxPlot([phi;phi],[ones(length(zeta1),1);2*ones(length(zeta2),1);3*ones(length(zeta3),1);4*ones(length(zeta4),1);5*ones(length(zeta),1)])
figure;notBoxPlot([zeta;zeta],[ones(length(zeta1),1);2*ones(length(zeta2),1);3*ones(length(zeta3),1);4*ones(length(zeta4),1);5*ones(length(zeta),1)])

figure;histogram(mu1,8,'BinLimits',[1000 5500])
figure;histogram(mu2,8,'BinLimits',[1000 5500])
figure;histogram(mu3,8,'BinLimits',[1000 5500])
figure;histogram(mu4,8,'BinLimits',[1000 5500])
figure;histogram(mu,8,'BinLimits',[1000 5500])

figure;histogram(phi1,8,'BinLimits',[0 2])
figure;histogram(phi2,8,'BinLimits',[0 2])
figure;histogram(phi3,8,'BinLimits',[0 2])
figure;histogram(phi4,8,'BinLimits',[0 2])
figure;histogram(phi,8,'BinLimits',[0 2])

figure;histogram(zeta1,8,'BinLimits',[0 2])
figure;histogram(zeta2,8,'BinLimits',[0 2])
figure;histogram(zeta3,8,'BinLimits',[0 2])
figure;histogram(zeta4,8,'BinLimits',[0 2])
figure;histogram(zeta,8,'BinLimits',[0 2])
