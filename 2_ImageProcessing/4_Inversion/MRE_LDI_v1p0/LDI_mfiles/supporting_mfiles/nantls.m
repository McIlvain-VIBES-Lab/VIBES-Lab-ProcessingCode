function X = nantls(A,B)
n = size(A,2); % n is the width of A (A is m by n)
C = [A B]; % C is A augmented with B.
nb = size(B,2);

if isempty(find(isnan(C))),
    [U S V] = svd(C,0); % find the SVD of C.
    VAB = V(1:n,1+n:end); % Take the block of V consisting of the first n rows and the n+1 to last column
    VBB = V(1+n:end,1+n:end); % Take the bottom-right block of V.
    X = -VAB/VBB;
else,
    X=NaN*ones(n,nb);
end;
