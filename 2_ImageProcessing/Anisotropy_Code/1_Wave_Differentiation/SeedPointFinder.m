function [SeedPoints] = SeedPointFinder(UT_mask,mask,allnvec)

[Surf_mask] = Surface_masker(mask);

List_all = [];
for ii = 1:size(mask,1)
    for jj = 1:size(mask,2)
        for kk = 1:size(mask,3)
            if Surf_mask(ii,jj,kk) == 1
                tmp = squeeze(UT_mask(ii,jj,kk,:));
                ix = find(tmp==max(tmp));

                for ll = 1:length(allnvec)
                    dotv(ll,1) = dot(allnvec(ix,:),allnvec(ll,:));
                end
                
                dtheta = acosd(dotv);
                
                [~,I] = sort(dtheta);
                dfx = cat(2,I,dtheta(I),tmp(I));
                ptx = find(dfx(:,2)>=50,1);
                tmpx = tmp(I);
                allnvecx = allnvec(I,:);
                
                wv = allnvecx(1:ptx,:).*repmat(tmpx(1:ptx),[1 3]);
                Dir = mean(wv,1)./norm(mean(wv,1));
                % A1 = norm(mean(wv,1));
                A = tmpx(1);
                x = ii;
                y = jj;
                z = kk;
                List_all = [List_all; [A, Dir, x, y, z]];
            else
            end
        end
    end
end

%1 Amplitude
%2-4 x, y, z vector components
%5-7 x, y, z location 

Orig_list = List_all;
SeedMap = zeros(size(mask));
SeedPoints = [];
trigger = 1;

while trigger == 1
    
Locat = find(List_all(:,1) == max(List_all(:,1)));
x = List_all(Locat,5);
y = List_all(Locat,6);
z = List_all(Locat,7);

SeedPoints = [SeedPoints; [List_all(Locat,5), List_all(Locat,6), List_all(Locat,7)]];
List_all(:,8) = zeros(size(List_all(:,1)));
List_all(:,9) = zeros(size(List_all(:,1)));

for ii = 1:length(List_all)
    List_all(ii,8) = acosd(dot([List_all(ii,2),List_all(ii,3),List_all(ii,4)],[List_all(Locat,2),List_all(Locat,3),List_all(Locat,4)]));
    List_all(ii,9) = sqrt((List_all(ii,5)-x)^2+(List_all(ii,6)-y)^2+(List_all(ii,7)-z)^2);
end

for ii = length(List_all):-1:1
    if List_all(ii,8) >= 45 && List_all(ii,9) >= 15
    else
        List_all(ii,:) = [];
    end
end


if length(List_all(:,1)) <= 200
    trigger = 0;
else
end

end

for ii = 1:size(mask,1)
    for jj = 1:size(mask,2)
        for kk = 1:size(mask,3)
            for ll = 1:size(SeedPoints,1)
                if SeedPoints(ll,1) == ii && SeedPoints(ll,2) == jj && SeedPoints(ll,3) == kk
                    SeedMap(ii,jj,kk) = 1;
                else
                end
            end
        end
    end
end

% figure;im(SeedMap)