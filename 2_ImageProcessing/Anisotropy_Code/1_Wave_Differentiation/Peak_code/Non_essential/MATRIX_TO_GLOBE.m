function [tmp] = MATRIX_TO_GLOBE(Matrix,x,y,z,N,radmat)
for z = 18:38
    Theta = radmat(:,1);
    Phi1 = radmat(:,2);
    tmp  = squeeze(Matrix(x,y,z,:,:));
    tmp_re = flip(reshape(flip(tmp,1),[1600 1]),1);
    Sphmat = [Theta Phi1 tmp_re];
    pointsize = 500;
    Phi = Phi1-pi()/2-pi()/N;
    
    % Phi_perm = permute(Phi,[2 1]);
    % Theta_perm = permute(Theta,[2 1]);
    % tmp_re_perm = permute(tmp_re,[2 1]);
    figure;PlotSphereIntensity(Theta,Phi,tmp_re);
    % figure;scatter(Sphmat(:,1),Sphmat(:,2),pointsize,Sphmat(:,3),'filled','square');caxis([0 .5]);
    
    % figure;sphere3d(tmp,0,2*pi(),-pi()/2,pi()/2,1,1);
end