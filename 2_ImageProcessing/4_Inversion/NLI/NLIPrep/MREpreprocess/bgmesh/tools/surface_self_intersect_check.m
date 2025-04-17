function r = surface_self_intersect_check(e,p)
% Checks if the triangular surface defined in 'e' and 'p' is
% self-intersecting
% r = 0 is it is not a self-intersecting surface otherwise r = 1

r = 0;

edges = [e(:,[1,2]);
       e(:,[1,3]);
       e(:,[2,3])];
edges = sort(edges,2);
edges = unique(edges,'rows');
facets_bbx = GetFacetsBBX(e,p);
% compute diagnoal of BBX of whole surface
gmin = min(p);
gmax = max(p);
d = norm(gmax-gmin);
mytol = d*1e-8;
for i=1:size(edges,1)
    n1 = edges(i,1);
    n2 = edges(i,2);
    p1 = p(n1,:);
    p2 = p(n2,:);
    
    [tf1 idx] = ismember(e,n1);
    tf1 = sum(tf1,2)>0;
    [tf2 idx] = ismember(e,n2);
    tf2 = sum(tf2,2)>0;
    tf = ~(tf1 & tf2);

    te = e;
    cmin = min(p1,p2);
    cmax = max(p1,p2);

    bf = facets_bbx(:,1) > cmax(1) & facets_bbx(:,2) > cmax(2) & ...
         facets_bbx(:,3) > cmax(3) & facets_bbx(:,4) < cmin(1) & ...
         facets_bbx(:,5) < cmin(2) & facets_bbx(:,6) < cmin(3);
    bf = ~bf & tf;
    te = te(bf,:);
    if isempty(te)
        continue
    end
    r = intersect_seg_triangles_mex(p1,p2,te,p,mytol);
    if r ~= 0
        break
    end
end
