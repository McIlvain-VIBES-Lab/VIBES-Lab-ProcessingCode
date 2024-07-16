function BW2 = centcluster(BW)
BW = circshift(BW,-5,2);
L = bwlabel(BW,4);
NumClusters = numel(unique(L)) -1;

Centers = zeros(NumClusters,2);
CenterLinIdices = zeros(NumClusters,1);

for k = 1:NumClusters

%// Find indices for elements forming each cluster.
    [r, c] = find(L==k);

%// Sort the elements to know hot many rows and columns the cluster is spanning.
    [~,y] = sort(r);
    c = c(y);
    r = r(y);

    NumRow = numel(unique(r));
    NumCol = numel(unique(c));

%// Calculate the approximate center of the cluster.
    CenterCoord = [r(1)+floor(NumRow/2) c(1)+floor(NumCol/2)];

%// Actually this array is not used here but you might want to keep it for future reference.
    Centers(k,:) = [CenterCoord(1) CenterCoord(2)];

%// Convert the subscripts indices to linear indices for easy reference.
if CenterCoord(1) > length(BW) 
    CenterCoord(1) = length(BW);
end
if CenterCoord(2) > length(BW)
    CenterCoord(2) = length(BW);
end
    CenterLinIdices(k) = sub2ind(size(BW),CenterCoord(1),CenterCoord(2));
end

%// Create output matrix full of 0s, except at the center of the clusters.
BW2 = false(size(BW));
BW2(CenterLinIdices) = 1;
BW2 = circshift(BW2,5,2);