function [x2 y2 z2 C2] = Sphericalise(x, y, z, C)

splitLine = 1; % set to 1 initially as we want it to do the loop at least
% once, and MATLAB doesn't have a do while loop

x2 = x(:);
y2 = y(:);
z2 = z(:);
C2 = C(:);

while splitLine
    T = convhulln([x2 y2 z2]);
    
    % It is best to split lines on a per line basis (and not a per triangle
    % basis) as each line is shared by two triangles, so any line that is to be
    % split will be split twice (and the new point will be in the same place)
    lines = unique([T(:,1) T(:,2); T(:,1) T(:,3); T(:,2) T(:,3)], 'rows');

    splitLine = 0;
    for line = lines'

        A = struct('x', x2(line(1)), 'y', y2(line(1)), 'z', z2(line(1)), 'C', C2(line(1)));
        B = struct('x', x2(line(2)), 'y', y2(line(2)), 'z', z2(line(2)), 'C', C2(line(2)));

        % if the length of the geodesic is longer than pi/10
        if acos(dot2(A, B)) > pi/10
            splitLine = 1;
            
            % add point between A and B
            
            % It is possible to calculate elevation and azimuth for A and
            % B, and then interpolate this to get elevation and azimuth for
            % P, and finally convert these back to cartesian coordinates,
            % but this 1) is inefficient and 2) doesn't work around the
            % poles (at least if theta and phi are linearly interpolated).

            P.x = (A.x + B.x) / 2;
            P.y = (A.y + B.y) / 2;
            P.z = (A.z + B.z) / 2;
            magnitude = sqrt(P.x^2 + P.y^2 + P.z^2);
            P.x = P.x / magnitude;
            P.y = P.y / magnitude;
            P.z = P.z / magnitude;
            
            P.C = (A.C + B.C) / 2;
            
            x2(end+1) = P.x;
            y2(end+1) = P.y;
            z2(end+1) = P.z;
            C2(end+1) = P.C;

        end
    end
end

end

function d = dot2(A, B)
d = dot([A.x A.y A.z], [B.x B.y B.z]);
end