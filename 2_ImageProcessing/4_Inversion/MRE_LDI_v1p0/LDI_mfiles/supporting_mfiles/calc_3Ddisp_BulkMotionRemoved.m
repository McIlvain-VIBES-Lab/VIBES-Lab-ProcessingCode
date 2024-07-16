function [Brfo_all,Bifo_all,X,Y,Z, coef,BulkMotion]=calc_3Ddisp_BulkMotionRemoved(Brf_all,Bif_all)
% written by RJO on 2/22/13
% revised 3/29/13 assuming all pixels outside ROI are NaN
% revised 11/18/14 to enforce all rigid body constraints 
% the coupled equations are
%  u_RB = u_R2*Y + u_R3*Z + u_R4
%  v_RB = -u_R2*X +v_R3*Z + v_R4
%  w_RB = -u_R3*X - v_R3*Y + w_R4
% matrix is assembled so that the 6 variables can be found using least
% squares
%  solution vector is {u_R2 u_R3 vR3 u_R4 v_R4 w_R4}
%  the first three values are dimensionless
%  the last three values are the rigid body translation w/ units of Brf_all
% processes a stack of displacement data together
if nargin < 2,
    disp('Processing real data');
else
    disp('Processing complex data');
end
 [X Y Z] = meshgrid(1:size(Brf_all,2), 1:size(Brf_all,1),1:size(Brf_all,3));
    newmat=zeros(3*length(X(:)),6);
    % fill in
    newmat(1:length(X(:)),1)=Y(:);
    newmat(1:length(X(:)),2)=Z(:);
    newmat(1:length(X(:)),4)=ones(length(X(:)),1);
    newmat(length(X(:))+(1:length(X(:))),1)= -X(:);
    newmat(length(X(:))+(1:length(X(:))),3)= Z(:);
    newmat(length(X(:))+(1:length(X(:))),5)= ones(length(Y(:)),1);
    newmat(2*length(X(:))+(1:length(X(:))),2)= -X(:);
    newmat(2*length(X(:))+(1:length(X(:))),3)= -Y(:);
    newmat(2*length(X(:))+(1:length(X(:))),6)= ones(length(Z(:)),1);  
    
    if nargin <2,
        newvec=Brf_all(:);
    else
        newvec=Brf_all(:)+1i*Bif_all(:);
    end
    userow=find(~isnan(newvec)>0);

    newmatuse=newmat(userow,:); newvecuse=newvec(userow);
    
    coef= newmatuse\newvecuse;
    upred = newmat*coef; 
    predmat=reshape(upred,size(Brf_all));
    Brfo_all = Brf_all- real(predmat);
    if nargin <2,
        Bifo_all=zeros(size(Brfo_all)); 
    else
        Bifo_all = Bif_all - imag(predmat);
    end    
    BulkMotion=predmat;
end