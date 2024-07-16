function [Gp, Gdp, C, nreTLS, nre, R2TLS, R, usepts] = matLSfit3DMS_v3d(u,v,w,del2u,del2v,del2w,rho,fms,Nxy,mask,delta,fracNAN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Gp, Gdp, C, nreTLS, nre, R2TLS, R, usepts] = matLSfit3DMS_v3d(u,v,w,del2u,del2v,del2w,rho,fms,Nxy,mask,delta,fracNAN)
%
% Inverting the reduced EOM (compressional wave negligable or removed prior)
%
% All data points within specified kernel (Nx,Ny,Nz) are fit to the coupled
% EOM in a least-squares sense. Nx must equal Ny in this program, hence
% Nxy.
%
% x = aG' - bG'' + Rc1
% y = bG' + aG'' + Ic2
%
% MULTISLICE VERSION OF matLSfit_v1a.m
% !! ALL SLICE DATA USED IN FIT !!
%
% INPUTS
% u,v,w     : complex displacement field (fundamental harmonic), u v w
%             components
% del2u,... : complex laplacian field (of fundatental harmonic), u v w
%             components
% mask   : mask (optional)
% delta  : 1/lambda (cf. Okamoto et al. 2011)
% fracNAN : fraction of pixels that can be NaN before the fit is not run
%          
%
% OUTPUTS:
% Gp      : Storage Modulus
% Gdp     : Loss Modulus
% C       : constants C(:,:,1) == Re, C(:,:,2) == Im
% R2      : correlation coeficient
% nre     : LS Normalized Residual Error
%
% Erik H. Clayton
% Washington University in St. Louis
% 11 October 2010; REV 15 October 2010; REV 12 January 2011-1D Multislice
% Rev by RJO 16 Feb 2011; Rev 21 Feb 2011 to use singular value
% decomposition
% 9 September 2011 - error on line 105: p(2) = -p(1)*my_xshift +my_yshift;
% Rev by RJO 18 Nov 2012 to fit 3-D displacement field
% Rev by RJO 20 Nov 2012 to rescale 3-D displacements by largest std dev
%                        before fitting and compute nre using OLS formula
% Rev by RJO 02 Aug 2013 to use NaN from center of fitting region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ELASTOcmap.mat;  %commented out on 11/11/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gp  = zeros(size(u,1),size(u,2));  Gdp = zeros(size(u,1),size(u,2));
C  = zeros(size(u,1),size(u,2),6);   
R = zeros(size(u,1),size(u,2)); R2TLS = zeros(size(u,1),size(u,2));
nre = zeros(size(u,1),size(u,2));  nreTLS = zeros(size(u,1),size(u,2)); 
usepts = zeros(size(u,1),size(u,2));
ww = fms*2*pi;

% input a mask is optional
if nargin <12,
    fracNAN=0.25;
elseif nargin <11,
    delta = 1;
end

if nargin >= 10,
    maskNAN = real(mask);
elseif nargin == 9,
    mask = ones(size(u));
    maskNAN = mask;
end
maskNAN(mask == 0) = NaN;

% for MULTISLICE ALLOWED
if size(u,4) == 1,
    fprintf('\n3D TLS INVERSION Multislice               : matLSfit3DMS_v3d.m');
    fprintf('\nKernel Fit, xy dim                          : %3.0f',2*Nxy+1);
    fprintf('\nKernel Fit, z  dim                          : %3.0f\n',size(u,3));
    [mmm,nnn,ppp] = size(u);%[m,n,~] = size(u);
    
    %u = real(u).*maskNAN +  imag(u).*maskNAN;
    temp = u.*maskNAN;
    
    for mm = 1+Nxy:mmm-Nxy,
        for nn = 1+Nxy:nnn-Nxy,
            
            if isnan(temp(mm,nn,(ppp-1)/2+1)) == 1  %%%%% Don't consider masked regions
                
            else
                % parsing data into chuncks then vectorizing
                % u component
                x  = real(u(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                y  = imag(u(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                a = real(del2u(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                b = imag(del2u(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                
                % v component
                m  = real(v(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                n  = imag(v(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                c = real(del2v(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                d = imag(del2v(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                
                % w component
                p  = real(w(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                q  = imag(w(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                e = real(del2w(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                f = imag(del2w(mm-Nxy:mm+Nxy, nn-Nxy:nn+Nxy,:));
                
                
                my_xvarU = [a(:)'+1i*b(:)'];
                my_xvarV = [c(:)'+1i*d(:)'];
                my_xvarW = [e(:)'+1i*f(:)'];
                my_yvarU = -rho*ww^2*[x(:)'+1i*y(:)'];
                my_yvarV = -rho*ww^2*[m(:)'+1i*n(:)'];
                my_yvarW = -rho*ww^2*[p(:)'+1i*q(:)'];
                % check if there are any NaN, discard point if specified fraction of points are NaN;
                xNANsum = sum(isnan(my_xvarU));
                yNANsum = sum(isnan(my_yvarU));
                if max(xNANsum,yNANsum) > fracNAN*length(my_yvarU),
                    pp = ones(1,4)*(NaN + NaN*1i);
                    r2TLStmp = NaN; nreTLStmp = NaN; 
                    use = NaN;
                else
                    use = (~isnan(my_xvarU) & ~isnan(my_yvarU) & ~isnan(my_xvarV) & ~isnan(my_yvarV) & ~isnan(my_xvarW) & ~isnan(my_yvarW)) == 1;
                    my_xshiftU = nanmean(my_xvarU(:));my_xsfU = nanstd(my_xvarU(:));
                    my_yshiftU = nanmean(my_yvarU(:));my_ysfU = nanstd(my_yvarU(:));
                    my_xshiftV= nanmean(my_xvarV(:));my_xsfV = nanstd(my_xvarV(:));
                    my_yshiftV = nanmean(my_yvarV(:));my_ysfV = nanstd(my_yvarV(:));
                    my_xshiftW = nanmean(my_xvarW(:));my_xsfW = nanstd(my_xvarW(:));
                    my_yshiftW = nanmean(my_yvarW(:));my_ysfW = nanstd(my_yvarW(:));
                    % center and scale data (all data must be scaled by
                    % same standard deviation, use largest of 3 disp
                    % directions
                    my_xsf=max([my_xsfU,my_xsfV,my_xsfW]);
                    my_ysf=max([my_ysfU,my_ysfV,my_ysfW]);
                    my_xfitU = (my_xvarU(use)-my_xshiftU)/my_xsf;
                    my_yfitU = (my_yvarU(use)-my_yshiftU)/my_ysf;
                    my_xfitV = (my_xvarV(use)-my_xshiftV)/my_xsf;
                    my_yfitV = (my_yvarV(use)-my_yshiftV)/my_ysf;
                    my_xfitW = (my_xvarW(use)-my_xshiftW)/my_xsf;
                    my_yfitW = (my_yvarW(use)-my_yshiftW)/my_ysf;

                    my_xfit = [my_xfitU(:);my_xfitV(:);my_xfitW(:)];
                    my_yfit = [my_yfitU(:);my_yfitV(:);my_yfitW(:)]; %make these col vectors
                    
                    % fit data using TOTAL least squares fn (tls.m)
                    pstar = delta*nantls(my_xfit,1/delta*my_yfit);
                    % calculate covariance of x and y
                    tmpcov = cov(my_xfit,my_yfit);
                    % calculate r^2 based on covariance matrix
                    r2TLStmp = (tmpcov(1,2)*tmpcov(2,1))/(tmpcov(1,1)*tmpcov(2,2));
                    nreTLStmp = sqrt(1-r2TLStmp);
                    % convert fitted polynomial back to unscaled values and get kernel fit estimates G' and G''
                    pp(1) = pstar(1)*my_ysf/my_xsf;
                    pp(2) = -pp(1)*my_xshiftU + my_yshiftU;
                    pp(3) = -pp(1)*my_xshiftV + my_yshiftV;
                    pp(4) = -pp(1)*my_xshiftW + my_yshiftW;
                end
                
                % separate constants
                Gp(mm,nn)  = real(pp(1));
                Gdp(mm,nn) = imag(pp(1));
                C(mm,nn,1) = real(pp(2));
                C(mm,nn,2) = imag(pp(2));
                C(mm,nn,3) = real(pp(3));
                C(mm,nn,4) = imag(pp(3));
                C(mm,nn,5) = real(pp(4));
                C(mm,nn,6) = imag(pp(4));
                
                % set R2 amd nre measurements for each region
                R2TLS(mm,nn) = r2TLStmp; nreTLS(mm,nn) = nreTLStmp; usepts(mm,nn) = sum(use);
 
                %NRE computation to match method used for OLS
 
                 endX = length(x(:));    endY = endX+length(y(:));
               % Lx = length(find(x(:))); Ly = length(find(y(:)));
                
                endM = endY+length(m(:)); endN = endM+length(n(:));
                %Lm = length(find(m(:)));  Ln = length(find(n(:)));
                
                endP = endN+length(p(:)); endQ = endP+length(q(:));
                %Lp = length(find(p(:)));  Lq = length(find(q(:)));
                
                % least-squares fit of entire anatomical region
                Ia = ones(size(a(:)));
                Za = zeros(size(a(:)));
                
                % all 8 free parameters (2 moduli plus constants)
                X1 = [a(:), -b(:),  Ia,Za,Za,Za,Za,Za;
                    b(:),  a(:), Za, Ia,Za,Za,Za,Za;
                    c(:), -d(:), Za,Za, Ia,Za,Za,Za;
                    d(:),  c(:), Za,Za,Za, Ia,Za,Za;
                    e(:), -f(:), Za,Za,Za,Za, Ia,Za;
                    f(:),  e(:), Za,Za,Za,Za,Za, Ia];

                 % r for each region
                M=reshape([real(pp(:))';imag(pp(:))'],8,1); %for consistency
                rsq(:,:,1) = corrcoef(-rho*ww^2*x(:), X1(1:endX,:)*M);
                rsq(:,:,2) = corrcoef(-rho*ww^2*y(:), X1(endX+1:endY,:)*M);
                rsq(:,:,3) = corrcoef(-rho*ww^2*m(:), X1(endY+1:endM,:)*M);
                rsq(:,:,4) = corrcoef(-rho*ww^2*n(:), X1(endM+1:endN,:)*M);
                rsq(:,:,5) = corrcoef(-rho*ww^2*p(:), X1(endN+1:endP,:)*M);
                rsq(:,:,6) = corrcoef(-rho*ww^2*q(:), X1(endP+1:endQ,:)*M);
                R(mm,nn) = mean(rsq(1,2,:),3);
                
               % max normalized residual error
                Xfit = X1(1:endX,:)*M; % real eqn fit
                Yfit = X1(endX+1:endY,:)*M; % imag eqn fit
                Mfit = X1(endY+1:endM,:)*M;
                Nfit = X1(endM+1:endN,:)*M;
                Pfit = X1(endN+1:endP,:)*M;
                Qfit = X1(endP+1:endQ,:)*M;
                               
                risX = (-rho*ww^2*x(:) - Xfit);
                risY = (-rho*ww^2*y(:) - Yfit);
                risM = (-rho*ww^2*m(:) - Mfit);
                risN = (-rho*ww^2*n(:) - Nfit);
                risP = (-rho*ww^2*p(:) - Pfit);
                risQ = (-rho*ww^2*q(:) - Qfit);
                
                num = [risX, risY, risM, risN, risP, risQ].^2;
                den = (rho*ww^2*[x(:), y(:), m(:), n(:), p(:), q(:)]).^2;
                
                % rjo addition 9 Sept 2011 (this will make nreTLS = nre)
                % den = (rho*ww^2*[x(:)-mean(x(:)), y(:)-mean(y(:)), m(:)-mean(m(:)), n(:)-mean(n(:)), p(:)-mean(p(:)), q(:)-mean(q(:))]).^2;
                                
                nre(mm,nn) = sqrt(sum(sum(num))/sum(sum(den)));
                % END NRE computation
            
            
            end %if 
            
            clear x y m n p q a b c d e f my_* use
        end
    end
else
     fprintf('\nONLY SINGLE SLICE OF DATA CAN BE INVERTED AT THIS TIME\n')
end
    return