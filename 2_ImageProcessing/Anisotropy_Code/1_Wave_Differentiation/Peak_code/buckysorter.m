load buckyinfo
orig = allnvec;
AA = zeros(size(allnvec));
nn = 1;
while isempty(orig) == 0
    movemat = find(orig(:,3) == max(orig(:,3))); %locate point(s) of highest z value
    move_mat = orig(movemat,:); % call all x,y,z of that point 
    if length(movemat)>1
        movemat2 = find(move_mat(:,2) == max(move_mat(:,2)));
        move_mat2 = orig(movemat2,:);
        AA(nn,:) = move_mat2;
        orig(movemat2,:) = [];
    else
        move_mat
        AA(nn,:) = move_mat;
        orig(movemat,:) = [];
    end
    nn = nn+1
end

uniq = unique(allnvec(:,3));

for NN = 1:length(uniq)
    movemat = find(orig(:,3) == uniq(NN));
    move_mat = orig(movemat,:);
    BB = move_mat
end