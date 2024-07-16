% load bucky300.mat
% 
% V1index = zeros(80,80,48);
% V2index = zeros(80,80,48);
% dot_buckyV1 = zeros(length(allnvec),1);
% dot_buckyV2 = zeros(length(allnvec),1);
% nn = 0;
% 
% for ix=1:80
%     for iy = 1:80
%         for iz = 1:48
%             if mask(ix,iy,iz)
%                 nn = nn+1;
%                 if rem(nn,500)
%                     nn
%                 end
%             for kk = 1:length(allnvec)
%                 dot_buckyV1(kk,1) = dot(squeeze(V1x(ix,iy,iz,:))',allnvec(kk,:));
%                 dot_buckyV2(kk,1) = dot(squeeze(V2x(ix,iy,iz,:))',allnvec(kk,:));
%             end
%              V1index(ix,iy,iz) = find(dot_buckyV1==min(dot_buckyV1),1);
%              V2index(ix,iy,iz) = find(dot_buckyV2==min(dot_buckyV2),1);
%             end
%         end
%     end
% end

load bucky300.mat
allnvec_mat = permute(repmat(allnvec,[1 1 80 80 48]),[3 4 5 2 1]);

V1index = zeros(80,80,48);
V2index = zeros(80,80,48);
dot_buckyV1 = zeros(80,80,48,300);
dot_buckyV2 = zeros(80,80,48,300);

for kk = 1:length(allnvec)
    dot_buckyV1(:,:,:,kk) = dot(V1x,squeeze(allnvec_mat(:,:,:,:,kk)),4);
    dot_buckyV2(:,:,:,kk) = dot(V2x,squeeze(allnvec_mat(:,:,:,:,kk)),4);
end

for ii = 1:80
    ii;
    for jj = 1:80
        for kk = 1:48
            V1index(ii,jj,kk) = find(abs(squeeze(dot_buckyV1(ii,jj,kk,:)))==min(abs(squeeze(dot_buckyV1(ii,jj,kk,:)))),1);
            V2index(ii,jj,kk) = find(abs(squeeze(dot_buckyV2(ii,jj,kk,:)))==min(abs(squeeze(dot_buckyV2(ii,jj,kk,:)))),1);
        end
    end
end