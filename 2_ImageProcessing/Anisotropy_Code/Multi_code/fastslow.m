function [m_s,m_f] = fastslow(A,xy,z,allnvec)
m_s = zeros(xy,xy,z,3,92);
m_f = zeros(xy,xy,z,3,92);
for nn = 1:length(allnvec(:,1))
    tic
    for i = 1:xy
        for j = 1:xy
            for k = 1:z
                m_s(i,j,k,:,nn) = cross(squeeze(allnvec(nn,:)),squeeze(A(i,j,k,:)));
                m_s(i,j,k,:,nn) = cross(squeeze(allnvec(nn,:)),squeeze(m_s(i,j,k,:,nn)));
            end
        end
    end
    toc
end