function [percent20,Unorms,perc20] = percent(A1x,A2x,U_s_V1,U_f_V1,U_s_V2,U_f_V2,mask)
Atot = A1x+A2x;
A1norm = A1x./Atot.*mask;
A2norm = A2x./Atot.*mask;

UV1 = U_s_V1+U_f_V1; UV2 = U_s_V2+U_f_V2;
U_s = U_s_V1+U_s_V2; U_f = U_f_V1+U_f_V2;

Us1norm = (U_s_V1./UV1).*A1norm; Uf1norm = (U_f_V1./UV1).*A1norm;
Us2norm = (U_s_V2./UV2).*A2norm; Uf2norm = (U_f_V2./UV2).*A2norm;
Unorms= cat(4,Us1norm,Uf1norm,Us2norm,Uf2norm);

Us1perc20 = (abs(Us1norm)>.2).*mask; Uf1perc20 = (abs(Uf1norm)>.2).*mask;
Us2perc20 = (abs(Us2norm)>.2).*mask; Uf2perc20 = (abs(Uf2norm)>.2).*mask;
perc20 = cat(4,Us1perc20,Uf1perc20,Us2perc20,Uf2perc20);

Us1perc20sum = sum(sum(sum(Us1perc20))); Uf1perc20sum = sum(sum(sum(Uf1perc20)));
Us2perc20sum = sum(sum(sum(Us2perc20))); Uf2perc20sum = sum(sum(sum(Uf2perc20)));
totalvox = sum(sum(sum(abs(U_s)>0)));
totalvoxf= sum(sum(sum(abs(U_f)>0)));

Us1percent20 = Us1perc20sum/totalvox; Uf1percent20 = Uf1perc20sum/totalvoxf;
Us2percent20 = Us2perc20sum/totalvox; Uf2percent20 = Uf2perc20sum/totalvoxf;

percent20 = [Us1percent20,Uf1percent20,Us2percent20,Uf2percent20];
