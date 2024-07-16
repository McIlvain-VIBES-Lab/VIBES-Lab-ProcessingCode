function [peak_data] = peakfinder4(matrix,N1,N2,maxpeaks)
cutoff = 1;
Mat = matrix;
Phi_sect = pi()/N2; Theta_sect = 2*pi()/N1;
maxmat = max(max(max(Mat)));
peaks = zeros(maxpeaks,1); peaks_loc = zeros(maxpeaks,2);
row_loc = zeros(maxpeaks,1); col_loc = zeros(size(row_loc));
onoff = 1; numpeaks = 0; n1 =0;
space = 3*floor(N1/20);

while onoff == 1
    n1 = n1+1;
    Matrix = Mat.*(Mat>=(maxmat.*cutoff));
    [matsub] = centcluster(Matrix);
    [matsub1] = centcluster(Matrix);
    if isempty(find(matsub))
    else
        tmp1 = size(matrix);
        n = 0; z = 0;
        tmprow = ones(tmp1(1),1);
        num_row = zeros(size(tmprow));  
        num_col = zeros(size(tmprow));
        for rw = 1:tmp1(1)
            tmprow(rw) = tmprow(rw) - isempty(find(matsub1(rw,:)));
            num_row(rw) = sum(matsub(rw,:));
            n = n+1;
        end
        row_loc = find(tmprow);
        for rl = 1:length(row_loc)
            tmpcol = ones(tmp1(2),1);
            row = row_loc(rl);
            if row == 1 || row == N2
                col_loc = N1/2;
            else
                for col = 1:tmp1(2)
                    z = z+1;
                    tmpcol(col) = tmpcol(col) - isempty(find(matsub1(row,col)));
                    num_col(col) = sum(matsub(rw,:));
                end
                col_loc = find(tmpcol);
            end
            for cl = 1:length(col_loc) 
                column = col_loc(cl);
                if rem(row,pi()/N2) == 0
                    row = row/Phi_sect;
                end
                %if length
                if row <= space
                    if column+space > N1 
                        if Mat(row,N1) ~= 0 && Mat(row,column-space) ~= 0 && Mat(row+space,column) ~= 0 && Mat(1,column) ~= 0 && Mat(row+space-1,column) && Mat(row+space-2,column) && Mat(row,1) ~= 0 && Mat(row,2) ~= 0
                            if row > 1+floor(N1/20)
                                numpeaks = numpeaks+1;
                                peaks(numpeaks) = Mat(row,column);
                                for jj = 1:row+space
                                    if column+space > N1
                                        for kk = column-space:N1
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    elseif column-space < 1
                                        for kk = 1:column+space
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    end
                                end
                                column = column*Theta_sect;
                                row = row*Phi_sect;
                                peaks_loc(numpeaks,:) = [column,row];
                            else
                                if sum(Mat(row+space,:) ~=0) == N1 && sum(Mat(row+space-1,:) ~=0) == N1 && sum(Mat(row+space-2,:) ~=0) == N1 && sum(Mat(row+space-3,:) ~=0) == N1
                                    numpeaks = numpeaks+1;
                                    peaks(numpeaks) = Mat(row,column);
                                    column = column*Theta_sect;
                                    row = row*Phi_sect;
                                    peaks_loc(numpeaks,:) = [column,row];
                                    matsub(1:space+1,:) = 1;
                                end
                            end
                        end
                    elseif column-space <= 0
                        if Mat(row,1) ~= 0 && Mat(row,column+space) ~= 0 && Mat(row+space,column) ~= 0 && Mat(1,column) ~= 0 && Mat(row+space-1,column) && Mat(row+space-2,column) && Mat(row,N1-1) ~= 0 && Mat(row,N1) ~= 0
                            if row > 1+floor(N1/20)
                                numpeaks = numpeaks+1;
                                peaks(numpeaks) = Mat(row,column);
                                for jj = 1:row+space
                                    if column+space > N1
                                        for kk = column-space:N1
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    elseif column-space < 1
                                        for kk = 1:column+space
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    end
                                end
                                column = column*Theta_sect;
                                row = row*Phi_sect;
                                peaks_loc(numpeaks,:) = [column,row];
                            else
                                if sum(Mat(row+space,:) ~=0) == N1 && sum(Mat(row+space-1,:) ~=0) == N1 && sum(Mat(row+space-2,:) ~=0) == N1 && sum(Mat(row+space-3,:) ~=0) == N1
                                    numpeaks = numpeaks+1;
                                    peaks(numpeaks) = Mat(row,column);
                                    column = column*Theta_sect;
                                    row = row*Phi_sect;
                                    peaks_loc(numpeaks,:) = [column,row];
                                    matsub(1:space+1,:) = 1;
                                end
                            end
                        end
                    elseif Mat(row,column+space) ~= 0 && Mat(row,column-space) ~= 0 && Mat(row+space,column) ~= 0 && Mat(1,column) ~= 0 && Mat(row,column+space-1) ~= 0 && Mat(row,column-space+1) ~= 0 && Mat(row+space-1,column) && Mat(row,column+space-2) ~= 0 && Mat(row,column-space+2) ~= 0 && Mat(row+space-2,column)
                        if sum(Mat(row+space,:) ~=0) == N1 && sum(Mat(row+space-1,:) ~=0) == N1 && sum(Mat(row+space-2,:) ~=0) == N1 && sum(Mat(row+space-3,:) ~=0) == N1
                            numpeaks = numpeaks+1;
                            peaks(numpeaks) = Mat(row,column);
                            column = column*Theta_sect;
                            row = row*Phi_sect;
                            peaks_loc(numpeaks,:) = [column,row];
                            matsub(1:space+1,:) = 1;
                        end
                    end
                elseif row >= N2+1-space
                    if column+space > N1 
                        if Mat(row,N1) ~= 0 && Mat(row,column-space) ~= 0 &&  Mat(row-space,column) ~= 0 && Mat(N2,column) ~= 0 && Mat(row-space+1,column) ~= 0 && Mat(row-space+2,column) ~= 0 && Mat(row,2) ~= 0 && Mat(row,1) ~= 0
                            if row < N2-floor(N2/20)
                                numpeaks = numpeaks+1;
                                peaks(numpeaks) = Mat(row,column);
                                for jj = row-space:N2
                                    if column+space > N1
                                        for kk = column-space:N1
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    elseif column-space < 1
                                        for kk = 1:column+space
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    end
                                end
                                column = column*Theta_sect;
                                row = row*Phi_sect;
                                peaks_loc(numpeaks,:) = [column,row];
                            else
                                if sum(Mat(row-space,:) ~=0) == N2 && sum(Mat(row-space+1,:) ~=0) == N2 && sum(Mat(row-space+2,:) ~=0) == N2 && sum(Mat(row-space+3,:) ~=0) == N2 
                                    numpeaks = numpeaks+1;
                                    matsub(N2+1-space:N2,:) = 1;
                                    peaks(numpeaks) = Mat(row,column);
                                    column = column*Theta_sect;
                                    row = row*Phi_sect;
                                    peaks_loc(numpeaks,:) = [column,row];
                                end
                            end
                        end
                    elseif column-space <= 0
                        if Mat(row,1) ~= 0 && Mat(row,column+space) ~= 0 && Mat(row-space,column) ~= 0 && Mat(N2,column) ~= 0 && Mat(row-space+1,column) ~=0 && Mat(row-space+2,column) ~=0 && Mat(row,N1-1) ~= 0 && Mat(row,N1) ~= 0
                            if row < N2-floor(N2/20)
                                numpeaks = numpeaks+1;
                                peaks(numpeaks) = Mat(row,column);
                                for jj = row-space:N2
                                    if column+space > N1
                                        for kk = column-space:N1
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    elseif column-space < 1
                                        for kk = 1:column+space
                                            [jj,kk];
                                            matsub(jj,kk) = 1;
                                        end
                                    end
                                end
                                column = column*Theta_sect;
                                row = row*Phi_sect;
                                peaks_loc(numpeaks,:) = [column,row];
                            else
                                if sum(Mat(row-space,:) ~=0) == N2 && sum(Mat(row-space+1,:) ~=0) == N2 && sum(Mat(row-space+2,:) ~=0) == N2 && sum(Mat(row-space+3,:) ~=0) == N2  
                                    numpeaks = numpeaks+1;
                                    matsub(N2+1-space:N2,:) = 1;
                                    peaks(numpeaks) = Mat(row,column);
                                    column = column*Theta_sect;
                                    row = row*Phi_sect;
                                    peaks_loc(numpeaks,:) = [column,row];
                                end
                            end
                        end
                    elseif Mat(row,column+space) ~= 0 && Mat(row,column-space) ~= 0 &&  Mat(row-space,column) ~= 0 && Mat(N2,column) ~= 0 && Mat(row,column+space-1) ~= 0 && Mat(row,column-space+1) ~= 0 &&  Mat(row-space+1,column) && Mat(row,column+space-2) ~= 0 && Mat(row,column-space+2) ~= 0 &&  Mat(row-space+2,column)
                        if sum(Mat(row-space,:) ~=0) == N2 && sum(Mat(row-space+1,:) ~=0) == N2 && sum(Mat(row-space+2,:) ~=0) == N2 && sum(Mat(row-space+3,:) ~=0) == N2  
                            numpeaks = numpeaks+1;
                            matsub(N2+1-space:N2,:) = 1;
                            peaks(numpeaks) = Mat(row,column);
                            column = column*Theta_sect;
                            row = row*Phi_sect;
                            peaks_loc(numpeaks,:) = [column,row];
                        end
                    end
                elseif column <= space
                    if Mat(row,column+space) ~= 0 && Mat(row+space,column) ~= 0 && Mat(row-space,column) ~= 0 && Mat(row,column+space-1) ~= 0 && Mat(row+space-1,column) ~= 0 && Mat(row-space+1,column) ~= 0 && Mat(row,column+space-2) ~= 0 && Mat(row+space-2,column) ~= 0 && Mat(row-space+2,column) ~= 0 
                        numpeaks = numpeaks+1;
                        overlap = space-column;
                        peaks(numpeaks) = Mat(row,column);
                        for jj = row-space:row+space
                            matsub(jj,1:column+space) = 1;
                            matsub(jj,N1-overlap) = 1;
                        end
                        column = column*Theta_sect;
                        row = row*Phi_sect;
                        peaks_loc(numpeaks,:) = [column,row];
                        
                    end
                elseif column >= N1-space
                    if Mat(row,column-space) ~= 0 && Mat(row+space,column) ~= 0 && Mat(row-space,column) ~= 0 && Mat(row,column-space+1) ~= 0 && Mat(row+space-1,column) ~= 0 && Mat(row-space+1,column) ~= 0 && Mat(row,column-space+2) ~= 0 && Mat(row+space-2,column) ~= 0 && Mat(row-space+2,column) ~= 0
                       numpeaks = numpeaks+1;
                       overlap = space+column-N1;
                       peaks(numpeaks) = Mat(row,column);
                       for jj = row-space:row+space
                           matsub(jj,column-space:N1) = 1;
                           matsub(jj,1:overlap) = 1;
                       end
                       column = column*Theta_sect;
                       row = row*Phi_sect;
                       peaks_loc(numpeaks,:) = [column,row];
                       
                    end
                else
                    if Mat(row,column+space) ~= 0 && Mat(row,column-space) ~= 0 && Mat(row+space,column) ~= 0 && Mat(row-space,column) ~= 0 && Mat(row,column+space-1) ~= 0 && Mat(row,column-space+1) ~= 0 && Mat(row+space-1,column) ~= 0 && Mat(row-space+1,column) ~= 0 && Mat(row,column+space-2) ~= 0 && Mat(row,column-space+2) ~= 0 && Mat(row+space-2,column) ~= 0 && Mat(row-space+2,column) ~= 0
                        numpeaks = numpeaks+1;
                        if column <= space
                            overlap = space-column;
                            for jj = row-space:row+space
                                matsub(jj,1:column+space) = 1;
                                matsub(jj,N1-overlap) = 1;
                            end
                        elseif column >= N1-space
                            overlap = space+column-N1;
                            for jj = row-space:row+space
                                matsub(jj,column-space:N1) = 1;
                                matsub(jj,1:overlap) = 1;
                            end
                        else
                            for jj = row-space:row+space
                                for kk = column-space:column+space
                                    [jj,kk];
                                    matsub(jj,kk) = 1;
                                end
                            end
                        end
                        peaks(numpeaks) = Mat(row,column);
                        column = column*Theta_sect;
                        row = row*Phi_sect;
                        peaks_loc(numpeaks,:) = [column,row];
                    end
                    if numpeaks == maxpeaks
                        onoff = 0;
                    end
                end
            end
            
        end
        Mat = Mat.*((Mat>0)-matsub);
    end   
    cutoff = cutoff-.001;
    if n1 == 1000
        break
    end
end

peak_data = cat(2,peaks_loc,peaks);