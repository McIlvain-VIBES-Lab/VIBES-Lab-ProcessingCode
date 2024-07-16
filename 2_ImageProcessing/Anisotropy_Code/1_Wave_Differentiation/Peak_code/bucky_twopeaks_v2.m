function [V1, V2, A1, A2] = bucky_twopeaks_v2(allnvec,data)

% load bucky200.mat

% points around peak = 7
% ptx = 42;
filtn = length(allnvec);

% peak 1
tmp = data; % (squeeze(UT_mask(37,26,28,:)));
ix = find(tmp==max(tmp));

for ii = 1:length(allnvec)
    dotv(ii,1) = dot(allnvec(ix,:),allnvec(ii,:));
end

dtheta = acosd(dotv);

[~,I] = sort(dtheta);
dfx = cat(2,I,dtheta(I),tmp(I));
ptx = find(dfx(:,2)>=50,1);

tmpx = tmp(I);
allnvecx = allnvec(I,:);

wv = allnvecx(1:ptx,:).*repmat(tmpx(1:ptx),[1 3]);
V1 = mean(wv,1)./norm(mean(wv,1));
% A1 = norm(mean(wv,1));
A1 = tmpx(1);

% peak 2
tmp = tmpx((ptx+1):end);
allnvec = allnvecx((ptx+1):end,:);
clear tmpx wv dfx npt dtheta dotv ix I Y

ix = find(tmp==max(tmp));

for ii = 1:length(allnvec)
    dotv(ii,1) = dot(allnvec(ix,:),allnvec(ii,:));
end

dtheta = acosd(dotv);

[~,I] = sort(dtheta);
dfx = cat(2,I,dtheta(I),tmp(I));
ptx = find(dfx(:,2)>=50,1);

tmpx = tmp(I);
allnvecx = allnvec(I,:);

wv = allnvecx(1:ptx,:).*repmat(tmpx(1:ptx),[1 3]);
V2 = mean(wv,1)./norm(mean(wv,1));
% A2 = norm(mean(wv,1));
A2 = tmpx(1);
%Add Comment Collapse



