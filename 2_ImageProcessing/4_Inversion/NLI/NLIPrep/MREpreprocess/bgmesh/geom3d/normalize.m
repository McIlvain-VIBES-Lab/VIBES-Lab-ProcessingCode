function n = normalize(v)
% function n = normalize(v)
% normalize the vector v, and return u.
%
% Songbai Ji (2/6/2006)
%
% Bug fixed when v is nXm matrix. (songbai 6/30/2006).

vv = v.*v;
vv = sqrt(sum(vv, 2));
n = v./repmat(vv,1,size(v,2));