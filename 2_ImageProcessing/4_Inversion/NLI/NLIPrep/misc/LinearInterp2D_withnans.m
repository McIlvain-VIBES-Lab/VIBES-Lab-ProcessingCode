function [ vout ] = LinearInterp2D_withnans(xin,yin,vin,xout,yout,cont)
%SplineInterp2D_withnans: Masked 2D cubic Linear interpolation - does not uses NaN values in the interpolation. 
%   Interpolates a 2D image using 1D Linear interpolation in
%   each of the 2 directions, one after the other. NaN values in the input
%   are not used in the interpolation, they are considered to be gaps
%   in the data. Interpolated values within these NaN sections are extrapolated 
%   from 'good' data and can be inaccurate. 
%   2 options are possible when there is more than section of good data
%   across a line in the image: 
%   1:  Continuous interpolation, where NaN's are considered bad
%       datapoints, but the end of one section and start of the next are
%       expected to be smoothly connected. Use contspline='continuous'
%   2:  Non-Continuous interpolation, where each section of data alond a line is treated
%       as independent, and indepenedent sections are fitted, the final 
%       solution is made up from a patchwork of these independent fits. This
%       is the default behavior.
%   The interpolated data is processed, setting any interpolated value which 
%   is closer to a NaN than a non-NaN datapoint to zero. A mask of 'safe' values
%   can easily be created from the output using safemask=~isnan(vout).
%
%   Requires the function LinearInterp1D_withnans
%
%   INPUTS
%    xin: X locations of input data. X directionis assumed to be the 1st
%         matlab dimension, i.e. moving from vin(1,1) to vin(2,1) is a
%         movement in the x direction.
%    yin: Y locations of input data. Y directionis assumed to be the 2nd
%         matlab dimension, i.e. moving from vin(1,1) to vin(1,2) is a
%         movement in the y direction.
%    vin: Input data, containing NaN's where the measurements are not
%         relevant. Size=[length(xin) length(yin)]
%    xout,yout: X and Y locations of interpolated data.
%    cont(optional): Set to 'continuous; to use a continuous fit
%          across a line of data with breaks in it. Otherwise, an indpenedent 
%          section is built for each continuous section of data.
%   OUTPUT
%    vout: Interpolated values of vin, at the points supplied by xout and yout.
%          Values of vin which are NaN's are not used in the interpolation.
%          If [xout(i),yout(j)] is closer to a NaN value in vin,
%          vout(i,j)=NaN. Note that this distance measurement is not an
%          absolute distance, because each dimension is processed
%          seperatly, it checks the the y distance first, and then the x distance. 
%
%   Author: Matt Mcgarry  -matthew.d.mcgarry@darmtouth.edu
%           Thayer School of Engineering
%           Dartmouth College, Hanover NH.
%     Date: 22 Dec 2011.
%
%   Example:
%   load Image.mat;
%   load Mask.mat
%   Im(mask~=1)=nan; % <- Im is just a masked 2D Image, values outside mask=NaN.
%   xin=1:size(Im,1);yin=1:size(Im,2);
%   xout=1:0.7:size(Im,1);yout=1:0.7:size(Im,2); % interpolate by factor of 0.7
%   %EXAMPLE 1: Process each continuous section along a line individually 
%   Im_interp=LinearInterp2D_withnans(xin,yin2,Im,xout,yout,'no')
%   %EXAMPLE 2: Process multiple continuous sections along a line as a continuous interpolation. 
%   Im_interp=LinearInterp3D_withnans(xin,yin,Im,xout,yout,'continuous')
%   safemask=~isnan(Im_interp).

%% Process inputs
if(nargin<6) % neither btyp or contspline is supplied
    cont='no';
end

if(numel(xin)>length(xin))
    error('xin must be a 1D vector')
end
if(numel(yin)>length(yin))
    error('yin must be a 1D vector')
end
if(numel(xout)>length(xout))
    error('xout must be a 1D vector')
end
if(numel(yout)>length(yout))
    error('yout must be a 1D vector')
end
    

% ensure vin is correct size
if(size(vin,1)~=length(xin))||(size(vin,2)~=length(yin))  
    error('vin must be size [length(xin) length(yin)]')
end

dim=size(vin);
dimout=[length(xout) length(yout)];

% x direction first
VoutY=zeros([dim(1) dimout(2)]);
for ii=1:dim(1)
    VoutY(ii,:)= LinearInterp1D_withnans(yin,vin(ii,:),yout,cont);
    VoutY(ii,:)= Checkfornans(yin,vin(ii,:),yout,VoutY(ii,:));        
end

% Now the y direction
vout=zeros([dimout(1) dimout(2)]);
for ii=1:dimout(1)
    vout(:,ii)=LinearInterp1D_withnans(xin,VoutY(:,ii),xout,cont);
    vout(:,ii)=Checkfornans(xin,VoutY(:,ii),xout,vout(:,ii));    
end
clear VoutY

end

function vout2= Checkfornans(xin,vin,xout,vout)
%Checkfornans: Checks 1D interpolated values for proximity to NaN in the
%input. If the interpolated value is closer to a NaN than real data, sets
%it to NaN.

vout2=vout;
for ii=1:length(xout)
    [y,i]=min(abs(xout(ii)-xin));
    if(isnan(vin(i)))
        vout2(ii)=nan;
    end
end

end

