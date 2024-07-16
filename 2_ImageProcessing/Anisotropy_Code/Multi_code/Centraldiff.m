function [StrainMatrix] = Centraldiff(xmotion,ymotion,zmotion,xy,z);
%always reload in the x,Y,Zmotion before running again or else you'll
%double the filter
%SGfilter
Xmotion_1 = sgolayfilt(xmotion,1,3);
Ymotion_1 = sgolayfilt(ymotion,1,3);
Zmotion_1 = sgolayfilt(zmotion,1,3);
%BPfilter
%Xmotion = 
for jj=1:xy    
    for ii = 1:xy
        for kk = 1:z
            %[ii jj kk]
            if jj == 1
                CDyx(jj,ii,kk) = (Xmotion_1(jj+1,ii,kk) - Xmotion_1(jj,ii,kk));
                CDyy(jj,ii,kk) = (Ymotion_1(jj+1,ii,kk) - Ymotion_1(jj,ii,kk));
                CDyz(jj,ii,kk) = (Zmotion_1(jj+1,ii,kk) - Zmotion_1(ii,ii,kk));               
            elseif jj == xy
                CDyx(jj,ii,kk) = (Xmotion_1(jj,ii,kk) - Xmotion_1(jj-1,ii,kk));
                CDyy(jj,ii,kk) = (Ymotion_1(jj,ii,kk) - Ymotion_1(jj-1,ii,kk));
                CDyz(jj,ii,kk) = (Zmotion_1(jj,ii,kk) - Zmotion_1(jj-1,ii,kk));         
            else
                CDyx(jj,ii,kk) = (Xmotion_1(jj+1,ii,kk) - Xmotion_1(jj-1,ii,kk))/2;
                CDyy(jj,ii,kk) = (Ymotion_1(jj+1,ii,kk) - Ymotion_1(jj-1,ii,kk))/2;
                CDyz(jj,ii,kk) = (Zmotion_1(jj+1,ii,kk) - Zmotion_1(jj-1,ii,kk))/2;
            end
        end
    end
end
for jj=1:xy    
    for ii = 1:xy
        for kk = 1:z
           if ii == 1
                CDxx(jj,ii,kk) = (Xmotion_1(jj,ii+1,kk) - Xmotion_1(jj,ii,kk));
                CDxy(jj,ii,kk) = (Ymotion_1(jj,ii+1,kk) - Ymotion_1(jj,ii,kk));
                CDxz(jj,ii,kk) = (Zmotion_1(jj,ii+1,kk) - Zmotion_1(jj,ii,kk));                
            elseif ii == xy
                CDxx(jj,ii,kk) = (Xmotion_1(jj,ii,kk) - Xmotion_1(jj,ii-1,kk));
                CDxy(jj,ii,kk) = (Ymotion_1(jj,ii,kk) - Ymotion_1(jj,ii-1,kk));
                CDxz(jj,ii,kk) = (Zmotion_1(jj,ii,kk) - Zmotion_1(jj,ii-1,kk));                
            else
                CDxx(jj,ii,kk) = (Xmotion_1(jj,ii+1,kk) - Xmotion_1(jj,ii-1,kk))/2;
                CDxy(jj,ii,kk) = (Ymotion_1(jj,ii+1,kk) - Ymotion_1(jj,ii-1,kk))/2;
                CDxz(jj,ii,kk) = (Zmotion_1(jj,ii+1,kk) - Zmotion_1(jj,ii-1,kk))/2; 
           end
        end
    end
end
for jj = 1:xy    
    for ii = 1:xy
        for kk = 1:z
           if kk == 1
                CDzx(jj,ii,kk) = (Xmotion_1(jj,ii,kk+1) - Xmotion_1(jj,ii,kk));
                CDzy(jj,ii,kk) = (Ymotion_1(jj,ii,kk+1) - Ymotion_1(jj,ii,kk));
                CDzz(jj,ii,kk) = (Zmotion_1(jj,ii,kk+1) - Zmotion_1(jj,ii,kk));               
            elseif kk == z
                CDzx(jj,ii,kk) = (Xmotion_1(jj,ii,kk) - Xmotion_1(jj,ii,kk-1));
                CDzy(jj,ii,kk) = (Ymotion_1(jj,ii,kk) - Ymotion_1(jj,ii,kk-1));
                CDzz(jj,ii,kk) = (Zmotion_1(jj,ii,kk) - Zmotion_1(jj,ii,kk-1));               
            else
                CDzx(jj,ii,kk) = (Xmotion_1(jj,ii,kk+1) - Xmotion_1(jj,ii,kk-1))/2;
                CDzy(jj,ii,kk) = (Ymotion_1(jj,ii,kk+1) - Ymotion_1(jj,ii,kk-1))/2;
                CDzz(jj,ii,kk) = (Zmotion_1(jj,ii,kk+1) - Zmotion_1(jj,ii,kk-1))/2;               
           end
        end
    end
end



StrainMatrix(:,:,:,1) = CDxx; StrainMatrix(:,:,:,2) = CDxy; StrainMatrix(:,:,:,3) = CDxz;
StrainMatrix(:,:,:,4) = CDyx; StrainMatrix(:,:,:,5) = CDyy; StrainMatrix(:,:,:,6) = CDyz; 
StrainMatrix(:,:,:,7) = CDzx; StrainMatrix(:,:,:,8) = CDzy; StrainMatrix(:,:,:,9) = CDzz;
%figure;im(StrainMatrix)