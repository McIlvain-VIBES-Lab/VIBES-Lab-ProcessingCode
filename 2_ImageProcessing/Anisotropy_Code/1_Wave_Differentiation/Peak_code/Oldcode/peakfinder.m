function [peaks_loc,num_peaks] = peakfinder(Sphmat,BW2)
matrix = flip(reshape(Sphmat(:,3),[50 50]),1);
peaks = matrix(BW2 == 1);
num_peaks = length(peaks);
peaks_loc = zeros(3,num_peaks);
for pp = 1:num_peaks
    if length(Sphmat(Sphmat(:,3)==peaks(pp))) == 50
        peaks_loc(:,pp) = [0 0 peaks(pp)];
    else
        tmp = Sphmat(cat(2,Sphmat(:,3)==peaks(pp),Sphmat(:,3)==peaks(pp),Sphmat(:,3)==peaks(pp)));
        peaks_loc(:,pp) = tmp(:,1);
    end
end
peaks_loc = permute(peaks_loc,[2 1]);
