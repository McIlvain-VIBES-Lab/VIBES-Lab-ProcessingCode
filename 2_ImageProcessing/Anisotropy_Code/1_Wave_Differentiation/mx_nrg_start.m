function [x,y,z] = mx_nrg_start(UT_mask,mask)

[Surf_mask] = Surface_masker(mask);
[Surf_mask2] = Surface_masker(Surf_mask);

Max_energy_all = max(UT_mask,[],4);
Surf_maxes = Max_energy_all.*Surf_mask2;
%tmp = max(max(max(Surf_maxes)));
locat = (Surf_maxes == max(max(max(tmp))));
for ii = 1:size(Surf_mask,1)
    for jj = 1:size(Surf_mask,2)
        for kk = 1:size(Surf_mask,3)
            if locat(ii,jj,kk) == 1
                x = ii;
                y = jj;
                z = kk;
            else
            end
        end
    end
end
