function plane = fit2Plane(X,Y,Z)
% function plane = fit2Plane(X,Y,Z);
% fit points into a plane, represented in parametric form
% Based on: %"Least Squares Fitting of Data" by David Eberly.
%
% Songbai Ji 6/25/2006.

X=X(:); Y=Y(:); Z=Z(:);
if (length(X) ~=length(Y) || length(X) ~= length(Z))
    error('length of X,Y, and Z must be same!');
end
if length(X) <3;
    error('# of points must not be less than 3!');
end

M = [sum(X.*X), sum(X.*Y), sum(X);
    sum(X.*Y), sum(Y.*Y), sum(Y);
    sum(X), sum(Y), length(X)];
N = [sum(X.*Z);
    sum(Y.*Z);
    sum(Z)];
ABC = M\N;

x = [0 1 1]; y = [0 0 1];
z = ABC(1)*x + ABC(2)*y + ABC(3);
plane = createPlane([x(1),y(1),z(1)], ...
    [x(2),y(2),z(2)], ...
    [x(3),y(3),z(3)]);