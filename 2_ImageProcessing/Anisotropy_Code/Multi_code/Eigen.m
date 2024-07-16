function [Eigvec,Eigval] = Eigen(Eps,xy,z)
for i =1:xy;
    for j = 1:xy;
        for k = 1:z;
            [i,j,k];
            eigtest(1,1) = Eps(i, j, k, 1, 1);
            eigtest(1,2) = Eps(i, j, k, 1, 2);
            eigtest(1,3) = Eps(i, j, k, 1, 3);
            eigtest(2,1) = Eps(i, j, k, 2, 1);
            eigtest(2,2) = Eps(i, j, k, 2, 2);
            eigtest(2,3) = Eps(i, j, k, 2, 3);
            eigtest(3,1) = Eps(i, j, k, 3, 1);
            eigtest(3,2) = Eps(i, j, k, 3, 2);
            eigtest(3,3) = Eps(i, j, k, 3, 3);
            [V, D] = eig(eigtest,'vector');
            tmp = D(1);
            tmp2 = V(:,1);
            Eigval(i,j,k) = tmp;
            Eigvec(i,j,k,:) = tmp2;
        end    
    end   
end
