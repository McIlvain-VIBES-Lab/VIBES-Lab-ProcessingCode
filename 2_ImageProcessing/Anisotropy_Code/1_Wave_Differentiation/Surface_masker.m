function [Surf_mask] = Surface_masker(mask)
Surf_mask = zeros(size(mask));

for ii = 2:size(mask,1)-1
    for jj = 2:size(mask,2)-1
        for kk = 2:size(mask,3)-1
             if mask(ii,jj,kk) == 1
                for mm = ii-1:ii+1
                    for nn = jj-1:jj+1
                        for oo = kk-1:kk+1
                            % mask(nn,mm,oo)
                            if mask(mm,nn,oo) == 0
                                Surf_mask(ii,jj,kk) = 1;
                            else
                            end
                        end
                    end
                end
            else
            end
        end
    end
end