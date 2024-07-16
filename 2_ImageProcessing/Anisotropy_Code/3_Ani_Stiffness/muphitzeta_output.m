[xa_AP1] = combo_2waves_LDI(1,1,1);
[xa_AP2] = combo_2waves_LDI(1,1,2);
[xa_both1] = combo_2waves_LDI(2,1,1);
[xa_both2] = combo_2waves_LDI(2,1,2);

[xc_AP1] = combo_2waves_LDI_crit(1,1,1);
[xc_AP2] = combo_2waves_LDI_crit(1,1,2);
[xc_both1] = combo_2waves_LDI_crit(2,1,1);
[xc_both2] = combo_2waves_LDI_crit(2,1,2);

save('combox.mat','xa_AP1','xa_AP2','xa_both1','xa_both2','xc_AP1','xc_AP2','xc_both1','xc_both2')
% save('combor2.mat','r2a_AP1','r2a_AP2','r2a_both1','r2a_both2','r2c_AP1','r2c_AP2','r2c_both1','r2c_both2')

mu_body = xa_AP1(:,1,1); mu_body1 = mu_body(mu_body > 0);
phi_body = xa_AP1(:,2,1)./mu_body; phi_body1 = phi_body(phi_body > 0);
zeta_body = xa_AP1(:,3,1)./mu_body; zeta_body1 = zeta_body(zeta_body > 0);
figure;histogram(mu_body1)
figure;histogram(phi_body1)
figure;histogram(zeta_body1)

mu_genu = xa_AP1(:,1,2); mu_genu1 = mu_genu(mu_genu > 0);
phi_genu = xa_AP1(:,2,2)./mu_genu; phi_genu1 = phi_genu(phi_genu > 0);
zeta_genu = xa_AP1(:,3,2)./mu_genu; zeta_genu1 = zeta_genu(zeta_genu> 0);

mu_spl = xa_AP1(:,1,3); mu_spl1 = mu_spl(mu_spl > 0);
phi_spl = xa_AP1(:,2,3)./mu_spl; phi_spl1 = phi_spl(phi_spl > 0);
zeta_spl = xa_AP1(:,3,3)./mu_spl; zeta_spl1 = zeta_spl(zeta_spl > 0);

mu_all = xa_AP1(:,1,4); mu_all1 = mu_all(mu_all > 0);
phi_all = xa_AP1(:,2,4)./mu_all; phi_all1 = phi_all(phi_all > 0);
zeta_all = xa_AP1(:,3,4)./mu_all; zeta_all1 = zeta_all(zeta_all > 0);