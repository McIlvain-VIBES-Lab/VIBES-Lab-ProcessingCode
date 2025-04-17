function dest = transformPoint3d(point, trans)
%TRANSFORMPOINT3D : tranform a point with a 3D affine transform
%
%   PT2 = TRANSFORMPOINT3D(PT1, TRANS).
%   where PT1 has the form [xp yp zp], and TRANS is a [3x3], [3x4], [4x4]
%   matrix, return the point transformed with affine transform TRANS.
%
%   Format of TRANS can be one of :
%   [a b c]   ,   [a b c j] , or [a b c j]
%   [d e f]       [d e f k]      [d e f k]
%   [g h i]       [g h i l]      [g h i l]
%                                [0 0 0 1]
%
%   PT2 = TRANSFORMPOINT3D(PT1, TRANS) also work when PTA is a [Nx3] array
%   of double. In this case, PT2 has the same size as PT1.
%
%   See also :
%   translation3d
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 10/02/2005.
%


% eventually add null translation
if size(trans, 2)==3
    trans = [trans zeros(size(trans, 1), 1)];
end

% eventually add normalization
if size(trans, 1)==3
    trans = [trans;0 0 0 1];
end

%old version, but bad result 
%res = [point ones(size(point, 1), 1)]*trans';

res = (trans*[point ones(size(point, 1), 1)]')';

dest = res(:,1:3);