function [nvec_new, theta_prime, Amp_new, output] = angle_compare(x, y, z, nvec, theta_prime_max, allnvec, UT_mask)

%load bucky200.mat
amplist = squeeze(UT_mask(x,y,z,:));

for ii = 1:length(allnvec)
    dotv(ii,1) = dot(nvec,allnvec(ii,:));
end

dtheta = acosd(dotv);

[~,I] = sort(dtheta); % Sort the dirctions based upon angle tothe 
dfx = cat(2,I,dtheta(I),amplist(I));

ptx = find(dfx(:,2)<=theta_prime_max,1,'last'); % Find last point in within the acceptable range of angles
% ptx_all = find(dfx(:,2)<=theta_prime_max); 

%figure;bar(dfx(:,2),dfx(:,3))
%figure;bar(dfx(ptx_all,2),dfx(ptx_all,3))

amplist_order = amplist(I); %Order the amplitude based upon the sorting
amplist_theta = amplist_order(1:ptx); %Truncate the sorted amplitudes
allnvec_order = allnvec(I,:); %Order the directions based upon the sorting
allnvec_theta = allnvec_order(1:ptx,:); %Truncate the sorted directions

wv = allnvec_theta.*repmat(amplist_theta,[1 3]);
nvec_new = mean(wv,1)./norm(mean(wv,1));
theta_prime = acosd(dot(nvec,nvec_new)); %Angle between input and new wave direction
Amp_fit = polyfit(dfx(1:ptx,2),amplist_theta,1);
Amp_new = Amp_fit(2); % Find Amplitude of new peak

if Amp_new <= .25 % If amplitude is too low aka noise, end calculation
    Amp_new = 0;
    nvec_new = [0,0,0];
    theta_prime = 0;
    output = -1;
else
    output = 2;
end
