function [ vout ] = LinearInterp3D_withnans(xin,yin,zin,vin,xout,yout,zout,cont)
%SplineInterp3D_withnans: Masked 3D cubic Linear interpolation - does not uses NaN values for the interpolation. 
%   Interpolates a 3D stack of data using 1D Linear interpolation in
%   each of the 3 directions, one after the other. NaN values in the input
%   are not used in the interpolation, they are considered to be gaps
%   in the data. Interpolated values within these NaN sections are extrapolated 
%   from 'good' data and often highly inaccurate. 
%   2 options are possible when there is more than section of good data
%   across a line: 
%   1:  Continuous interpolation, where NaN's are considered bad
%       datapoints, but the end of one section and start of the next are
%       expected to be smoothly connected. Use cont='continuous'
%   2:  Non-Continuous interpolation, where each section of data along a line is treated
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
%         matlab dimension, i.e. moving from vin(1,1,1) to vin(2,1,1) is a
%         movement in the x direction.
%    yin: Y locations of input data. Y directionis assumed to be the 2nd
%         matlab dimension, i.e. moving from vin(1,1,1) to vin(1,2,1) is a
%         movement in the y direction.
%    zin: Z locations of input data. Z directionis assumed to be the 3rd
%         matlab dimension, i.e. moving from vin(1,1,1) to vin(1,1,2) is a
%         movement in the z direction.
%    vin: Input data, containing NaN's where the measurements are not
%         relevant. Size=[length(xin) length(yin) length(zin)]
%    xout,yout,zout: X, Y and Z locations of interpolated data.
%    cont(optional): Set to 'continuous; to use a continuous interpolation
%          across a line of data with breaks in it. Otherwise, an indpenedent 
%          section is built for each continuous section of data.
%   OUTPUT
%    vout: Interpolated values of vin, at the points supplied by xout,yout,zout.
%          Values of vin which are NaN's are not used in the interpolation.
%          If [xout(i),yout(j),zout(k)] is closer to a NaN value in vin,
%          vout(i,j,k)=NaN. Note that this distance measurement is not an
%          absolute 3D distance, because each dimension is processed
%          seperatly, it checks the z distance first, then the y distance,
%          then the x distance. 
%
%   Author: Matt Mcgarry  -matthew.d.mcgarry@darmtouth.edu
%           Thayer School of Engineering
%           Dartmouth College, Hanover NH.
%     Date: 22 Dec 2011.
%           Update 6/5/2017: stopped extrapolating past the edges of the 
%           input data.  
%
%   Example:
%   load MRE_3DMotionData.mat;
%   load Mask.mat
%   Ur=A(:,:,:,1).*cos(P(:,:,:,1));
%   Ur(mask~=1)=nan; % <- Ur is just a masked 3D stack of data, values outside mask=NaN.
%   xin=1:size(Ur,1);yin=1:size(Ur,2);zin=1:size(Ur,3);
%   xout=1:0.7:size(Ur,1);yout=1:0.7:size(Ur,2);zout=1:0.7:size(Ur,3);
%   %EXAMPLE 1: Process each continuous section along a line individually 
%   Ur_interp=LinearInterp3D_withnans(xin,yin,zin,Ur,xout,yout,zout,'no')
%   %EXAMPLE 2: Process multiple continuous sections along a line as a
%   continuous section
%   Ur_interp=LinearInterp3D_withnans(xin,yin,zin,Ur,xout,yout,zout,'continuous')
%   safemask=~isnan(Ur_interp).

%% Process inputs
if(nargin<8) % neither btyp or contspline is supplied
    cont='no';
end

if(numel(xin)>length(xin))
    error('xin must be a 1D vector')
end
if(numel(yin)>length(yin))
    error('yin must be a 1D vector')
end
if(numel(zin)>length(zin))
    error('zin must be a 1D vector')
end
if(numel(xout)>length(xout))
    error('xout must be a 1D vector')
end
if(numel(yout)>length(yout))
    error('yout must be a 1D vector')
end
if(numel(zout)>length(zout))
    error('zout must be a 1D vector')
end
    

% ensure vin is correct size
if(size(vin,1)~=length(xin))||(size(vin,2)~=length(yin))||(size(vin,3)~=length(zin))
    error('vin must be size [length(xin) length(yin) length(zin)]')
end

dim=size(vin);
dimout=[length(xout) length(yout) length(zout)];

% Slice direction first
VoutZ=zeros([dim(1) dim(2) dimout(3)]);
for ii=1:dim(1)
    for jj=1:dim(2)
        VoutZ(ii,jj,:)= LinearInterp1D_withnans(zin,vin(ii,jj,:),zout,cont);
        VoutZ(ii,jj,:)= Checkfornans(zin,vin(ii,jj,:),zout,VoutZ(ii,jj,:));
        
    end
end

% Now the second index
VoutYZ=zeros([dim(1) dimout(2) dimout(3)]);
for ii=1:dim(1)
    for jj=1:dimout(3)
        VoutYZ(ii,:,jj)=LinearInterp1D_withnans(yin,VoutZ(ii,:,jj),yout,cont);
        VoutYZ(ii,:,jj)=Checkfornans(yin,VoutZ(ii,:,jj),yout,VoutYZ(ii,:,jj));
    end
end
clear VoutZ

% Now the final index
vout=zeros([dimout(1) dimout(2) dimout(3)]);
for ii=1:dimout(2)
    for jj=1:dimout(3)
        vout(:,ii,jj)=LinearInterp1D_withnans(xin,VoutYZ(:,ii,jj),xout,cont);
        vout(:,ii,jj)=Checkfornans(xin,VoutYZ(:,ii,jj),xout,vout(:,ii,jj));
    end
end
clear VoutYZ


% Set all values outside of the data bounds to NaN
vout(xout<min(xin),:,:)=0;
vout(:,yout<min(yin),:)=0;
vout(:,:,zout<min(zin))=0;
vout(xout>max(xin),:,:)=0;
vout(:,yout>max(yin),:)=0;
vout(:,:,zout>max(zin))=0;



end

function vout2= Checkfornans(xin,vin,xout,vout)
%Checkfornans: Checks 1D interpolated values for proximity to NaN in the
%input. If the interpolated value is closer to a NaN than real data, sets
%it to NaN. Also sets to nan if it is too close to the edges. 

vout2=vout;
for ii=1:length(xout)
    [y,i]=min(abs(xout(ii)-xin));
    if(isnan(vin(i)))
        vout2(ii)=nan;
    end
end

end

