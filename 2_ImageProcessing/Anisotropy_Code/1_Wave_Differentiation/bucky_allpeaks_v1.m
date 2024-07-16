function [pointsall, points, A, V] = bucky_allpeaks_v1(UT_mask,allnvec,x,y,z)
% load bucky200.mat
jj = 0;
trig = 1;
% points around peak = 7
% ptx = 42;
tmp = squeeze(UT_mask(x,y,z,:));

while trig == 1
    jj = jj+1;
clear tmpx wv dfx npt dtheta dotv ix I Y
filtn = length(allnvec);
% peak 1
ix = find(tmp==max(tmp));

for ii = 1:length(allnvec)
    dotv(ii,1) = dot(allnvec(ix,:),allnvec(ii,:));
end

dtheta = acosd(dotv);

[~,I] = sort(dtheta);
dfx = cat(2,I,dtheta(I),tmp(I));
ptx = find(dfx(:,2)>=50,1);
figure;bar(dfx(:,2),dfx(:,3))
tmpx = tmp(I);
allnvecx = allnvec(I,:);

wv = allnvecx(1:ptx,:).*repmat(tmpx(1:ptx),[1 3]);
allnvecs(jj,:) = allnvec(ix,:);
V(jj,:) = mean(wv,1)./norm(mean(wv,1));
% A1 = norm(mean(wv,1));
A(jj) = tmpx(1);

% peak 2
tmp = tmpx((ptx+1):end);
allnvec = allnvecx((ptx+1):end,:);

if (max(A)/A(jj))>4
    trig = 0;
    A(jj) = [];
    V(jj,:) = [];
    allnvecs(jj,:) = [];
else
    trig = 1;
end
end

load bucky200.mat
points = zeros(length(allnvec(:,1)),1);
for ii = 1:length(allnvecs)
    points = points+(allnvec(:,1)==allnvecs(ii,1));
end
pointsall = find(points);

disp(sprintf('%i peaks found for point %i,%i,%i',ii,x,y,z))
