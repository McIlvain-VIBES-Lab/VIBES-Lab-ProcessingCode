function [peaks_loc,num_peaks] = peakfinder2(radmat,Rnew,BW2,N)
matrix = flip(reshape(Rnew,[N N]),1);
peaks = matrix(BW2 == 1);
num_peaks = length(peaks);
peaks_loc = zeros(2,num_peaks);
for pp = 1:num_peaks
    if length(radmat(Rnew==peaks(pp))) > 1
        peaks_loc(:,pp) = [0 0];
    else
        tmp = radmat(cat(2,Rnew==peaks(pp),Rnew==peaks(pp)));
        peaks_loc(:,pp) = tmp(:,1);
    end
end
peaks_loc = permute(peaks_loc,[2 1]);
