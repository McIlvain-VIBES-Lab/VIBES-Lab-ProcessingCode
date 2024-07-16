function [peaks,peaks_loc,peak_data] = peakfinder3(matrix,cutoff,peaksnum,N)
Phi_sect = pi()/N;
Theta_sect = 2*pi()/N;
space = floor(N/10);
peaks = zeros(peaksnum,1);
row_loc = zeros(peaksnum,1);
col_loc = zeros(size(row_loc));
Mat = matrix;
for ii = 1:peaksnum
    Matrix = Mat.*(Mat>(max(max(max(Mat)))).*cutoff); %Matrix only showing points on the matrix over a cutoff percentage of the max
    % Matrx(i1,j1,k1,:,:) = Matrix;
    matsub = zeros(size(matrix));
    [matsub] = centcluster(Matrix); % find center of peaks over the cutoff
    peaks(ii) = Matrix(matsub == 1); % value of the peaks
    tmp1 = size(matrix);
    n = 0; z = 0;
    tmprow = ones(tmp1(1),1);
    tmpcol = ones(tmp1(2),1);
    for row = 1:tmp1(1)
        tmprow(row) = tmprow(row) - isempty(find(matsub(row,:)));
        n = n+1;
    end
    row_loc(ii) = find(tmprow);
    for col = 1:tmp1(2)
        z = z+1;
        tmpcol(col) = tmpcol(col) - isempty(find(matsub(:,col)));
    end
    col_loc(ii) = find(tmpcol);
    if row_loc(ii) == N
        %1
        matsub(N-space:N,:) = 1;
    elseif row_loc(ii) == 1
        matsub(1:space+1,:) = 1;
        %2
    elseif col_loc(ii) <= space
        %3
        overlap = col_loc(ii) - space;
        for jj = row_loc(ii)-space:row_loc(ii)+space
            matsub(jj,1:col_loc(ii)+space) = 1;
            matsub(jj,N+1-overlap) = 1;
        end
    elseif col_loc(ii) > N-space
        %4
        overlap = space+col_loc(ii)-N;
        for jj = row_loc(ii)-space:row_loc(ii)+space
            matsub(jj,col_loc(ii)-space:N) = 1;
            matsub(jj,1:overlap) = 1;
        end
    else
        %5
        for jj = row_loc(ii)-space:row_loc(ii)+space
            for kk = col_loc(ii)-space:col_loc(ii)+space
                [jj,kk]
                matsub(jj,kk) = 1;
            end
        end
    end
    Mat = Mat.*(ones(size(Mat))-matsub);
end

col_loc = col_loc*Theta_sect;
row_loc = row_loc*Phi_sect;
peaks_loc = [col_loc,row_loc];
peak_data = [col_loc,row_loc,peaks];

%% peaks_loc = permute(peaks_loc,[2 1]);
