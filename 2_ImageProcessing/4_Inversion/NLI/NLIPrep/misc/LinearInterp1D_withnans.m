function [yout]=LinearInterp1D_withnans(xin,yin,xout,cont)
%LinearInterp1D_withnans: Perfroms a 1D Linear interpolation, does
%not use NaN values. The function can fit different sections of data 
% independently, the first continuous section is found, and fitted. The 
% remaining section is then recursively passed back into the function to 
% independently fit each section. Note that data outside of the sections is
% VERY unreliable, more than about 1 data spacing away. If there is only 1
% datapoint in a section, a constant value is returned for yout in the
% vicinity of that section.
% To fit a single function for all non-NaN data, use contspline='continuous'.
% There is potential for errors when there is multiple data sections because 
% the two are linked by a linear funtion, whereas they may be independent in 
% reality, so end effects may affect interpolation around the edge of breaks 
% in the data.
% INPUTS:
%   xin: x coordinates of the data - 1D vector, any orientation
%   yin: y values of the data, containing nan's where data is not present. - 1D vector, any orientation
%   xout: x cordinates of the desired interpolated values. - 1D vector, any orientation
%   cont(optional): set to 'continuous' to fit a single function to all
%           non-NaN datapoints. Otherwise seperate sections are fitted for
%           each data section 
%
%   Author: Matt Mcgarry  -matthew.d.mcgarry@dartmouth.edu
%           Thayer School of Engineering
%           Dartmouth College, Hanover NH.
%     Date: 22 Dec 2011.
%
% Examples: 
% xin=1:0.8:25;
% yin=sin(xin);
% xout=1:0.1:25;
% yin([1:2 12:16 end-2:end])=nan; % Create some gaps in the data
%
% EXAMPLE 1: independent interpolation of each section
% yout=LinearInterp1D_withnans(xin,yin,xout,'no');
% figure
% plot(xin,yin,'ro:',xout,yout,'b.')
%
%
% EXAMPLE 2: natural spline boundary conditions, continuous interpolation
% yout=LinearInterp1D_withnans(xin',yin,xout',1,'continuous');
% figure
% plot(xin,yin,'ro:',xout,yout,'b.')
% Note that the interpolation is continuous across regions of NaN's in the 
% data. 


%% Process inputs
if(nargin<4) % cont is not supplied
    cont='no';
end

% ensure output is same shape as xout;
yout=zeros(size(xout));

%Requires row vectors 
if(numel(xin)>length(xin))
    error('xin must be a 1D vector')
else
    xin=reshape(xin,[1 length(xin)]);
end

if(numel(yin)>length(yin))
    error('yin must be a 1D vector')
else
    yin=reshape(yin,[1 length(yin)]);
end

if(length(xin)~=length(yin))
    error('xin and yin must be same length')
end


%% Start the interpolation


if(sum(isnan(yin))==numel(yin)) % all entries are NaNs, return NaN.
    yout(:)=nan;
elseif(sum(~isnan(yin))==1) % If there is only one datapoint, return yout as constant
    yout=ones(size(yout))*yin(~isnan(yin));    
else % Process as normal 
    


if(strcmp(cont,'continuous'))
    Yspl=yin(~isnan(yin));
    Xspl=xin(~isnan(yin));
    yout=interp1(Xspl,Yspl,xout,'linear','extrap');
else % fit independent functions for each section of data
    
    % Account for up to 2 seperate sections here. - recursive call sorts out
    % other sections.
    hasvals=~isnan(yin);
    ind=1:length(yin);
    lim=[min(ind(hasvals)) max(ind(hasvals))];
    
    if(sum(isnan(yin(lim(1):lim(2))))==0) % only 1 section, treat same as continuous.
        Yspl=yin(~isnan(yin));
        Xspl=xin(~isnan(yin));
        yout=interp1(Xspl,Yspl,xout,'linear','extrap');
    else % More than one section, Fit indepenedent functions
        % Find end point of first section and start of second
        firstnan=lim(1);
        nextval=lim(1);
        foundnan=false;
        foundnext=false;
        ii=lim(1);
        while((ii<=lim(2))&&(~foundnan||~foundnext))
            if(isnan(yin(ii))&&~foundnan)
                firstnan=ii;
                foundnan=true;
            end

            if(~isnan(yin(ii))&&foundnan)
                nextval=ii;
                foundnext=true;
            end
            ii=ii+1;
        end
       
         % Fit spline to First section
        Yspl=yin(lim(1):firstnan-1);
        Xspl=xin(lim(1):firstnan-1);
        if(length(Yspl)==1) % catch case where there is only 1 datapoint in first section - return constant answer.
            Yint1=ones(size(yout))*Yspl(1);
        else            
            Yint1=interp1(Xspl,Yspl,xout,'linear','extrap');            
        end
        
        
        % Recursivly call function to fit independent functions to second
        % section, which may be be composed of other sections
        yinsec2=yin;
        yinsec2(1:firstnan)=nan;        
        [Yint2]=LinearInterp1D_withnans(xin,yinsec2,xout,cont);
                
        Xend1=xin(firstnan-1);
        Xstart2=xin(nextval);
        % Assign Uint value from whichever Section is closer to the end of
        % section 1 or start of section 2.
        for ii=1:length(Yint1)
            if abs(xout(ii)-Xend1)<abs(xout(ii)-Xstart2)
                yout(ii)=Yint1(ii);
            else
                yout(ii)=Yint2(ii);
            end
        end  
    end   
end

end

%% Figure to Check Results
if(0==1) % Show Figures for testing
    figure
    plot(xin,yin,'ro');
    hold on
    plot(xout,yout,'b.-')
    ylim([1.1*min(yin) 1.1*max(yin)])
end


end

