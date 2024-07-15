 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y
%   Brad Sutton, Univ. Michigan, June 2002

if a.is.empty
	error empty
end

G = a.G;
we = a.we;

i = sqrt(-1);
tt = a.tt;

L = size(a.int,1);

if ~(size(tt) == size(a.int,2))
  sprintf('Number of time points not the same between P and Int')
  keyboard
end

if (L == 1)
   tau = 0;
else
  tau = (max(tt)-min(tt)+eps)/(L-1);
end
Int = a.int;
TE = min(tt);
sizeG = size(G);


if ~a.is.transpose 
        vi = exp(-i*we(:)*TE).*vi(:);
        vo_tmp = zeros(sizeG(1),L);                 
         BB = G;
        for ll = 1:L                                 
            Wo = exp(-i*col(we)*((ll-1)*tau));
            vo_tmp(:,ll) = BB*(Wo.*col(vi));
        end              
        Int = permute(Int, [2 1]);
        aa = Int;
        vo = sum(vo_tmp.*aa, 2);       
        if a.flgswt
            vo = vo.*a.swt;
        end	           

else        
        if a.flgswt
	    vi = vi.*a.swt;  % Transpose of real sinc-weighting
        end
        vo_tmp = zeros(sizeG(2), L);
        BB = G;
        for ll = 1:L
            Wo = exp(i*conj(we)*((ll-1)*tau));
            aa_per_shot = Int(ll,:)';
            aa = aa_per_shot;
            vo_tmp(:,ll) = Wo(:).*(BB'*(aa.*col(vi)));              
        end
        vo = sum(vo_tmp, 2);
        vo = exp(i*conj(we(:))*TE).*vo;
	    %vo = vo(:);
end





