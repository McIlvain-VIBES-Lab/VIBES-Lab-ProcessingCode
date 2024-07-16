addpath('300filt')
load BuckyOut.mat
load bucky300.mat
n1 = zeros(80,80,48,3);
n2 = zeros(80,80,48,3);
Ufilt1 = zeros(80,80,48,3);
Ufilt2 = zeros(80,80,48,3);
m_s1x = zeros(80,80,48,3);
m_f1x = zeros(80,80,48,3);
m_s2x = zeros(80,80,48,3);
m_f2x = zeros(80,80,48,3);

for ii = 1:80
    for jj = 1:80
        for kk = 1:48
            tic
            filt1 = V1index(ii,jj,kk);
            filt2 = V2index(ii,jj,kk);
            n1 = squeeze(allnvec(filt1,:)');
            n2 = squeeze(allnvec(filt2,:)'); 
            A = squeeze(dtiall(ii,jj,kk,:));
            tmp1 = load(sprintf('U_filt%d.mat',filt1)); U_filt1 = tmp1.U_filt;
            tmp2 = load(sprintf('U_filt%d.mat',filt2)); U_filt2 = tmp2.U_filt;
            Ufilt1(ii,jj,kk,:) = squeeze(U_filt1(ii,jj,kk,:));
            Ufilt2(ii,jj,kk,:) = squeeze(U_filt2(ii,jj,kk,:));
            m_s1 = cross(n1,A)/norm(cross(n1,A));
            m_f1 = cross(n1,m_s1)/norm(cross(n1,m_s1));
            m_s2 = cross(n2,A)/norm(cross(n2,A));
            m_f2 = cross(n2,m_s2)/norm(cross(n2,m_s2));
            m_s1x(ii,jj,kk,:) = m_s1;
            m_f1x(ii,jj,kk,:) = m_f1;
            m_s2x(ii,jj,kk,:) = m_s2;
            m_f2x(ii,jj,kk,:) = m_f2;
            toc
        end
    end
end

% m_s1x = repmat(permute(m_s1,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);
% m_s2x = repmat(permute(m_s2,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);
% m_f1x = repmat(permute(m_f1,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);
% m_f2x = repmat(permute(m_f2,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);

Us1 = repmat(dot(Ufilt1,m_s1x,4),[1 1 1 3]).*m_s1x; Uf1 = repmat(dot(Ufilt1,m_f1,4),[1 1 1 3]).*m_f1x; 
Us2 = repmat(dot(Ufilt2,m_s2x,4),[1 1 1 3]).*m_s2x; Uf2 = repmat(dot(Ufilt2,m_f2,4),[1 1 1 3]).*m_f2x; 

