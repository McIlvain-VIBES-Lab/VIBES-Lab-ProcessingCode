function [V1, V2, A1, A2] = bucky_twopeaks(data,X,cutoff)

%load bucky300.mat
allnvec = X;

% points around peak = 7
%ptx = 42;
filtn = length(allnvec);

% peak 1
tmp = data; % (squeeze(UT_mask(37,26,28,:)));
ix = find(tmp==max(tmp));

for ii = 1:filtn
    dotv(ii,1) = dot(allnvec(ix,:),allnvec(ii,:));
end

dtheta = acosd(dotv);

[~,I] = sort(dtheta);
dfx = cat(2,I,dtheta(I),tmp(I));

tmpx = tmp(I);
allnvecx = allnvec(I,:);
n = 0;
k = 0;
dfx2 = [0,0,0];
dfx3 = [0,0,0];

for ii = 1:length(dfx)
    if dfx(ii,2) >= cutoff
        n = n+1;
        dfx2(n,:) = dfx(ii,:);
    end
    if dfx(ii,2) >= cutoff-15
        k = k+1;
        dfx3(k,:) = dfx(ii,:);
    end
end

aver = mean(dfx2(:,3),1);
stddev = std(dfx2(:,3),1);
n = 0; 
for ii = 1:length(dfx2)
    if dfx2(ii,2) >= aver+stddev
        n = n+1;
        dfx_cut(n,:) = dfx2(ii,:);
    end
end

wv = allnvecx(1:ptx,:).*repmat(tmpx(1:ptx),[1 3]);
V1 = mean(wv,1)./norm(mean(wv,1));
% A1 = norm(mean(wv,1));
A1 = tmpx(1);

% peak 2
tmp = tmpx((ptx+1):end);
allnvec = allnvecx((ptx+1):end,:);
clear tmpx wv dfx npt dtheta dotv ix I Y

ix = find(tmp==max(tmp));

for ii = 1:(filtn-ptx)
    dotv(ii,1) = dot(allnvec(ix,:),allnvec(ii,:));
end

dtheta = acosd(dotv);

[~,I] = sort(dtheta);
dfx = cat(2,I,dtheta(I),tmp(I));

tmpx = tmp(I);
allnvecx = allnvec(I,:);

wv = allnvecx(1:ptx,:).*repmat(tmpx(1:ptx),[1 3]);
V2 = mean(wv,1)./norm(mean(wv,1));
% A2 = norm(mean(wv,1));
A2 = tmpx(1);