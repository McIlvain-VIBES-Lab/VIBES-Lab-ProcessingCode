function [h]=CompareConvergence(flist,truth)

nf=length(flist);

for ifile=1:nf
    load(flist(ifile).name);
    f(ilist).conv=conv;
    f(ilist).reconind=reconind;    
end

if(nargin<1)
    load

