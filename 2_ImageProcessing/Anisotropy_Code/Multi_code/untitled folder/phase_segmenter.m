function [X,Y,Z] = phase_segmenter(Xmotion,Ymotion,Zmotion);
for t = 1:8
    n = (t-1)*pi/4;
    Xcom = Xmotion*exp(i*n);
    X(:,:,:,t) = real(Xcom);
    Ycom = Ymotion*exp(i*n);
    Y(:,:,:,t) = real(Ycom);
    Zcom = Zmotion*exp(i*n);
    Z(:,:,:,t) = real(Zcom);
end
