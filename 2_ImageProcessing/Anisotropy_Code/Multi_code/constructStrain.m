function [Eps] = constructStrain(DiffMatrix)

Eps(:,:,:,1,1) = DiffMatrix(:,:,:,1);
Eps(:,:,:,1,2) = .5 * (DiffMatrix(:,:,:,2)+DiffMatrix(:,:,:,4));
Eps(:,:,:,1,3) = .5 * (DiffMatrix(:,:,:,3)+DiffMatrix(:,:,:,7));
Eps(:,:,:,2,1) = .5 * (DiffMatrix(:,:,:,4)+DiffMatrix(:,:,:,2));
Eps(:,:,:,2,2) = DiffMatrix(:,:,:,5);
Eps(:,:,:,2,3) = .5 * (DiffMatrix(:,:,:,6)+DiffMatrix(:,:,:,8));
Eps(:,:,:,3,1) = .5 * (DiffMatrix(:,:,:,7)+DiffMatrix(:,:,:,3));
Eps(:,:,:,3,2) = .5 * (DiffMatrix(:,:,:,8)+DiffMatrix(:,:,:,6));
Eps(:,:,:,3,3) = DiffMatrix(:,:,:,9);




