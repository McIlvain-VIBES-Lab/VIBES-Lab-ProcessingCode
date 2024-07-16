function [Matrx,radmat,cartmat,Ux,Uy,Uz,T0,P0] = spheregrid(Bifo_all,Brfo_all,mask)
tic
% U3abs = abs(U3); EnerU3 = sqrt(U3abs(:,:,:,1).^2+U3abs(:,:,:,2).^2+U3abs(:,:,:,3).^2);
N = 40;
Phi_sect = pi()/N; Theta_sect = 2*pi()/N;
R = ones(N^2,1);
Theta_numbers = ones(1,N);
Phi_numbers = -N/2;
for ii = 2:N
    Theta_numbers = cat(2,Theta_numbers,repmat(ii,1,N));
end
for jj = Phi_numbers+1:(N/2-1)
    Phi_numbers = cat(2,Phi_numbers,jj);
end
size(Phi_numbers);
Phi_numbers = repmat(Phi_numbers,1,N);
Phi = permute(Phi_sect*Phi_numbers,[2 1]);
Theta = permute(Theta_sect*Theta_numbers,[2 1]);
Phi1 = Phi+pi()/2+pi()/N;
radmat = [Theta Phi1]; sphmat = [Theta Phi1 R];
[X,Y,Z] = sph2cart(Theta,Phi,R);
cartmat = [X Y Z];


toc
[Ux,Uy,Uz,T0,P0,U3filt] = buckyanalyzer_CLJfilter(Bifo_all,Brfo_all,mask,cartmat);
pointsize = 500;

U3abs = abs(U3filt); EnerU3 = sqrt(U3abs(:,:,:,1).^2+U3abs(:,:,:,2).^2+U3abs(:,:,:,3).^2);

% Multx = squeeze(Ux(47,ii,28,:));
% Multy = squeeze(Uy(47,ii,28,:));
% Multz = squeeze(Uz(47,ii,28,:));
% Mult_tot = cat(2,abs(Multx),abs(Multy),abs(Multz));
% Mult = mean(Mult_tot,2);
% Rnew = Mult.*R;
% Sphmat = cat(2,radmat,Rnew);
% 
% figure;scatter(Sphmat(:,1),Sphmat(:,2),pointsize,Sphmat(:,3),'filled','square');caxis([0 1])
% figure;PlotSphereIntensity(rad2deg(abs(Sphmat(:,1))),rad2deg(abs(Sphmat(:,2))),rad2deg(abs(Sphmat(:,3))))
% 
% matrix = flip(reshape(Rnew,[50 50]),1);
% Matrix = matrix.*(matrix>.35);
% 
% [BW2] = centcluster(Matrix);
% [peaks_loc,num_peaks] = peakfinder(Sphmat,BW2);

peaks = zeros(80,80,48,10,2);
peak_number = zeros(80,80,48);
Matrx = zeros(80,80,48,N,N);
Mat_full = zeros(size(Matrx));
n = 0;
tic
for i1 = 1:80
    for j1 = 1:80
        for k1 = 1:48
            if mask(i1,j1,k1) == 1
                n = n+1;
                if rem(n,100) == 0
                    n
                    toc
                end
                %[i1 j1 k1]
                Multx = squeeze(Ux(i1,j1,k1,:));
                Multy = squeeze(Uy(i1,j1,k1,:));
                Multz = squeeze(Uz(i1,j1,k1,:));
                Mult_tot = cat(2,abs(Multx),abs(Multy),abs(Multz));
                Mult = sqrt(sum(Mult_tot.^2,2)); % EnerU3_2500
                Mult_norm = Mult/EnerU3(i1,j1,k1);
                %toc
                % Rnew = Mult.*R;
                Rnew = Mult_norm.*R;
                % Sphmat = cat(2,radmat,Rnew);
                matrix = flip(reshape(Rnew,[N N]),1);
                Mat_full(i1,j1,k1,:,:) = matrix;
                % Matrix = matrix.*(matrix>(max(max(max(matrix)))*.9));
                Matrix = matrix.*(matrix>0.2);
                [BW2] = centcluster(Matrix);
                Matrx(i1,j1,k1,:,:) = Matrix;
                % [peaks_loc,num_peaks] = peakfinder(Sphmat,BW2);
                [peaks_loc,num_peaks] = peakfinder2(radmat,Rnew,BW2,N);
                for peak = 1:num_peaks
                    peaks(i1,j1,k1,peak,:) = peaks_loc(peak,:);
                end
                peak_number(i1,j1,k1) = num_peaks;
                %toc
            end
        end
    end
end
toc

save('Peaks_Output.mat','peaks','peak_number','Matrx','Mat_full','-v7.3')


tic
save('Sphereoutput.mat','Sphmat','Mult','Rnew')
toc
% 
% cd .. 
% cd ..
% load('dtiall.mat')
% cd Tracts
% tract = load_nii('tractsCR.nii'); TractCR = tract.img>0; dti = (1-(dtiall(:,:,:,3)==0));
% CR = TractCR.*dti;
% figure;im(CR)
