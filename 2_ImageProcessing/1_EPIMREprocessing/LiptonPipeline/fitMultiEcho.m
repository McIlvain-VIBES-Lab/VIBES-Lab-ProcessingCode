% Author: Roman Fleysher
% Albert Einstein College of Medicine
% Columbia University Medical Center
%
% module load MATLAB/R2022b
%-----------------------------------------------
function fitMultiEcho(outDecayNiiFile, outFieldMapNiiFile, teColumnFile, inNiiMagnitude, inNiiPhase, inSeedRegion)
% fit R2 and field map if phase data is provided
% inNiiPhase is optional input (in radians)
% inSeedRegion is mask where we believe field is within bandwidth and carrierOffset will be set to 0
%              if inSeedRegion is non-zero everywhere, carrier offset is never applied
%              pixels where inSeedRegion =0 are processed iteratively
%              with carrier offset taken from neighboring pixel where it was previously determined
% inNiiPhase and inSeedRegion must be supplied as a pair
%
% output is always single float
% output outDecayNiiFile NIFTI file with these images:
%   optimal T2 contrast (zero in pixels where fit is impossible)
%   fitted spin density, can be NaN. Not anymore?
%   decay rate constant (in supplied time units), can be NaN. Not anymore? 
%   inverse error, if multiplied by decay rate and divided by noise in individual input magnitude, will be SNR on the decay rate
%
% output outFieldMapNiiFile (produced only if inNiiPhase, inSeedRegion are supplied !!!) NIFTI file with these images:
%   optimal T2 contrast (zero in pixels where fit is impossible)
%   fitted spin density, can be NaN. Not anymore?
%   field map (in supplied time units), if phase is supplied 
%   inverse error, if noise in individual input magnitude is divided by this, one gets error on field map in supplied time units
%   this inverse error (and validity of field map) could be different in future implementation from that in decay rate
%   for example, if input has two disjoint regions, then heterodyne carrier offset can not be carried over the gap of no magnitude data
%   Disjoint regions must be processed separately, by bringing information how frequency is offset between them


  % these reduce memory usage a little
  clearvars -global, close all
  % limit number of cores to 1
  % it shortens run time and stops enormous CPU usages
  maxNumCompThreads(1);

  assert(nargin == 4 || nargin == 6 , 'Supply filenames: outDecay, outFieldMap, teColumnFile, inMagnitudeTimeseries, optional pair: inPhaseTimeseries, inSeedRegion');

  if (~isdeployed)
    oldPath = path;
    addpath('~/NIFTI');
  end

  % load list of echo times from a text file
  fileID = fopen(teColumnFile,'r');
  assert(fileID ~= -1, 'Can not open file with TE list');
  teList = fscanf(fileID,'%f');
  fclose(fileID);
  fileID=[];

  % fit decay and get auxiliary variables (although they do not seem to speed up field map calculation)
  % echo times can be in arbitrary order for decay fit
  % for fitting field map, will have to order by TE because phase is limited to +/-2pi
  [decayNii, intensityScale, denominator, totalWeight, weightedX ] = internal_fitMultiEchoMagnitude( teList, inNiiMagnitude);
  save_untouch_nii(decayNii, outDecayNiiFile);

  if ( nargin == 6 )
    % perform iterative calculation of field map using carrier offset
    % Start with the seed region (maybe a voxel) in the middle of the brain where we think the frequency is definitely within the bandwidth.
    % Fit field there using carrier offset =0 (standard fit). Use fitted field value as carrier offset to fit neighboring pixels,
    % thus expanding region where field is known. On the next pass, use newly fitted values on the edge of the fitted region
    % as pixel-wise carrier offset to expand once more. Keep doing this until all pixels are fit.
    % Since there are several neighbors, take carrier offset from one of the neighbors with the largest fit weight
    % (large amplitude, slowest decaying and therefore more trust-worthy neighbor because noise in all pixels is the same).

    % this code is inefficient in two ways: it keeps reading input files on every iteration (to reduce memory usage)
    % and it should fit only edge pixels on each iteration instead of fitting everywhere and saving only the edge.
    % Need to pass edge map to internal_fitMultiEcho


    % load seed region.
    % should be 0 everywhere and 1 where carrier offset =0 is OK.
    % no anatomy
    seedRegion = load_untouch_nii(inSeedRegion);
    seedRegion.img = seedRegion.hdr.dime.scl_slope * single(seedRegion.img) + seedRegion.hdr.dime.scl_inter;

    seedRegion.img = seedRegion.img(:,:,:,1); 
    seedRegion.img(seedRegion.img ~=0 ) = 1;
    validField = seedRegion.img;

    assert(sum(size(totalWeight) == size(validField)) == 3, 'inSeedRegion and inMagnitudeTimeseries must have the same size');

    % if seed region is too small, neighboring pixels will have less than 6 to choose carrier offset from
    % if seed region is too small (less that 3x3x3 pixels) expand it
    % this is not correct solution because we do not check pixel are contiguous.
    % this has to be done outside? 
    if ( sum(validField(:)) < 6)
      for x=[-1:1]
      for y=[-1:1]
      for z=[-1:1]
        validField = validField + circshift(seedRegion.img, [x,y,z]);
      end
      end
      end

      validField(validField ~=0 ) = 1;

    end


    % fit using carrier offset=0 everywhere
    fieldMap = internal_fitMultiEchoPhase(teList, inNiiMagnitude, inNiiPhase, zeros(size(validField)), intensityScale, denominator, totalWeight, weightedX );
    fitWeight = decayNii.img(:,:,:,4);

    % loop until all pixels become valid
    % not all pixels might have valid data (anatomy is burned out or brain mask)
    % need to loop until valid region stops expanding?
    % or count pixels in anatomy?
    % input data does not have anatomy
    % while ( sum(validField(:)) < prod(size(carrierOffset)) )

    % loop until valid region stops expanding
    oldValidFieldCount = 0;
    while ( oldValidFieldCount < sum(validField(:) > 0))

      oldValidFieldCount = sum(validField(:));
      % set fitWeight to 0 outside of valid region
      fitWeight(validField == 0 ) = 0;
      % set to 0 at edge of FOV so that we do not wrap around carrier offset
      fitWeight(1,:,:)=0;
      fitWeight(:,1,:)=0;
      fitWeight(:,:,1)=0;
      fitWeight(end,:,:)=0;
      fitWeight(:,end,:)=0;
      fitWeight(:,:,end)=0;

      % use field map (over seed region) as carrier offset
      carrierOffset = fieldMap;
      carrierOffset(validField == 0) = 0;
      bestWeight = fitWeight;
      choiceCounter = validField ;

      % move by +/-1 pixel in all directions to find best neighbor
      for x=[-1:1]
      for y=[-1:1]
      for z=[-1:1]

        % choiceCounter will count how many neighbors were used to select carrier
        % minimum required is 6 for 3D data or 3 for 2D data
        % this is why seed region can not be too small
        choiceCounter = choiceCounter + circshift(validField, [x,y,z]);

        shiftedFitWeight     = circshift(fitWeight, [x,y,z]);
        shiftedCarrierOffset = circshift(carrierOffset, [x,y,z]);

        carrierOffset(shiftedFitWeight > bestWeight) = shiftedCarrierOffset(shiftedFitWeight > bestWeight);
        bestWeight(shiftedFitWeight > bestWeight)    = shiftedFitWeight(shiftedFitWeight > bestWeight);

      end
      end
      end

      % minimum required is 6 for 3D data or 3 for 2D data
      % if seed region is too small, these cuts will not enlarge valid pixels and while will loop forever
      validField(bestWeight > 0) = 1;
      validField(choiceCounter < 6)  = 0;

      fieldMap = internal_fitMultiEchoPhase(teList, inNiiMagnitude, inNiiPhase, carrierOffset, intensityScale, denominator, totalWeight, weightedX);
      fitWeight = decayNii.img(:,:,:,4);
    end % while

    %nii.img(:,:,:,2) = validField;
    %nii.img(:,:,:,2) = bestWeight;
    %nii.img(:,:,:,3) = carrierOffset;

    % apply validField to output anatomy
    decayNii.img(validField==0) = 0;
    decayNii.img(:,:,:,3) = fieldMap;
    %decayNii.img(:,:,:,4) = decayNii.img(:,:,:,4);

    decayNii.hdr.dime.dim(5) = 4;

    save_untouch_nii(decayNii, outFieldMapNiiFile);

  end % if ( nargin == 6 )

  if (~isdeployed)
    path(oldPath);
  end

end
%-----------------------------------------------

function [nii, intensityScale, denominator, totalWeight, weightedX ] = internal_fitMultiEchoMagnitude(teList, inNiiMagnitude)
  % fit T2 
  % output is always single float
  % output NIFTI file with these images:
  %     image with optimal T2 contrast (zero in pixels where fit is impossible)
  %     fitted spin density, should not be NaN
  %     decay rate, should not be NaN (in supplied time units)
  %     inverse error in decay rate, if multiplied by rate and divided by noise in individual magnitude, will be SNR in decay rate. Small is bad
  % intensityScale   used to divide pixel values to bring all within (0,1] range
  %                  without it, weight (given by square if intensity) could be too large
  %                  inside this function used only to rescale weights
  %                  the following must be multiplied by square of it
  % denominator      of both decay and error estimate (is 0 if decay can not be determined)
  % totalWeight      total weight, w, sum of magnitude^2
  % weightedX        sum of w*TE

  assert(nargin == 2, 'Supply filenames: teColumnFile, inMagnitudeTimeseries');

  %-------------------------------
  % load magnitude data
  % We will load one volume at a time to conserve memory.
  % get matrix size from header
  %-------------------------------
  fullNiiHeader = load_nii_hdr(inNiiMagnitude); %nii.hdr

  nX  = fullNiiHeader.dime.dim(2);
  nY  = fullNiiHeader.dime.dim(3);
  nZ  = fullNiiHeader.dime.dim(4);
  nT  = fullNiiHeader.dime.dim(5);

  assert(nT >= 2, 'At least 2 time points on input are expected');

  %-------------------------------
  % fit decay into logarithm of data
  % model: S=A*exp(r*t)
  %        y=slope*x + intercept
  %  slope = r, intercept = log(A),
  %  y=log(S), x=t
  %
  % assume noise in S is independent of time
  % therefore, noise in y is 1/S and weight in least squares is S^2
  % weighted least squares:
  % Y = H * beta with weights W = diag (1/sigma^2)
  %          +    -1  +
  % beta = (H W H)  (H W Y)
  % covariance:
  %            +    +     -1
  % beta * beta = (H W H )
  %
  % in our (linear) case
  %
  % intercept =  ( sum (w*x^2) * sum (w*y)   - sum (w*x) * sum (w*x*y) ) / denominator
  % rate      =  ( sum (w)     * sum (w*x*y) - sum (w*x) * sum (w*y)   ) / denominator
  % where
  % denominator = sum (w) * sum (w*x^2) - (sum (w*x))^2
  %-------------------------------

  assert(length(teList) == nT, 'Number of TEs differs from number of time points in NIFTI file');

  totalWeight = zeros(nX,nY,nZ); % sum (w)
  weightedX   = zeros(nX,nY,nZ); % sum (w*x) (x is TE)
  weightedX2  = zeros(nX,nY,nZ); % sum (w*x^2)

  weightedY   = zeros(nX,nY,nZ); % sum (w*y)
  weightedXY  = zeros(nX,nY,nZ); % sum (w*x*y)

  % as we loop, adjust intensityScale to be applied to weights only !
  for k = 1:nT

    % load one time point at a time
    nii = load_untouch_nii(inNiiMagnitude, k);

    % Since we use load_untouch_nii, no scale/intercept is
    % applied. Apply them and convert to floating number to allow NaN
    nii.img = nii.hdr.dime.scl_slope * single(nii.img) + nii.hdr.dime.scl_inter;

    % add small offset to avoid log(0). Can find better solution?
    % y = log(nii.img +0.00001);
    assert(min(nii.img(:)) >= 0, 'Input magnitude data can not be negative');

    % adjust intensityScale for weights
    if ( k == 1 )
      intensityScale = max(nii.img(:));
    end
    scaleFactor = max(nii.img(:)) / intensityScale;
    if (scaleFactor > 1)
      totalWeight = totalWeight / scaleFactor / scaleFactor;
      weightedX   = weightedX   / scaleFactor / scaleFactor;   
      weightedX2  = weightedX2  / scaleFactor / scaleFactor;   

      weightedY   = weightedY   / scaleFactor / scaleFactor;
      weightedXY  = weightedXY  / scaleFactor / scaleFactor;

      intensityScale = intensityScale * scaleFactor;
    end

    % replace NaN and +/-Inf in y with zero because added weight of that time point will be zero
    % NaN in decay should only be created if less than 2 time points have non-zero data.
    y = log(nii.img);
    y(isnan(y)) = 0;
    y(isinf(y)) = 0;

    weight = (nii.img / intensityScale) .* (nii.img / intensityScale);

    % elements of weighted least squares
    totalWeight = totalWeight + weight;
    weightedX   = weightedX   + weight*teList(k);
    weightedX2  = weightedX2  + weight*teList(k)*teList(k);   

    weightedY   = weightedY   + weight.*y;
    weightedXY  = weightedXY  + weight.*teList(k).*y;

  end;


  %-------------------------------
  % compute decay and 
  % optimal T2 contrast for output anatomy
  %-------------------------------
  denominator = totalWeight .* weightedX2 - weightedX .* weightedX;

  decayRate   = (totalWeight.*weightedXY - weightedX.*weightedY) ./ denominator; % decay per time unit
  intercept   = (weightedY.*weightedX2   - weightedX.*weightedXY) ./ denominator; % log(spin density)
  spinDensity = exp(intercept);
  inverseFitError  = intensityScale * sqrt(denominator ./ totalWeight);

  % if totalWeight == 0, then there is no data, set to zeros
  % if denominator <= 0, then only one time point useful (there is no data), set to zeros
  % denominator could be negative due to round-off errors
  notEnoughData = (denominator <= 0) | (totalWeight == 0);
  decayRate(notEnoughData) = 0;
  spinDensity(notEnoughData) = 0;
  inverseFitError(notEnoughData) = 0;

  optimalContrast = zeros(size(spinDensity));
  weightNorm      = zeros(size(spinDensity));

  for k = 1:nT

    nii = load_untouch_nii(inNiiMagnitude, k);
    nii.img = nii.hdr.dime.scl_slope * single(nii.img) + nii.hdr.dime.scl_inter;

    weight = teList(k).*spinDensity.*exp(decayRate*teList(k)) / intensityScale;
    optimalContrast = optimalContrast + nii.img.*weight;
    weightNorm = weightNorm + weight.*weight;

  end;

  optimalContrast(notEnoughData ~=1) = optimalContrast(notEnoughData ~=1) ./ sqrt(weightNorm(notEnoughData ~=1));

  % this is not needed
  % replace NaN or Inf with zero to indicate where data not available
  % this should happen only if data from less than two distinct TEs is available
  % optimalContrast(isnan(optimalContrast) | isinf(optimalContrast)) = 0;
  % optimalContrast(notEnoughData) = 0;

  %-------------------------------
  % output anatomy, spin density and decay time constant
  %-------------------------------
  % force save as float
  nii.hdr.dime.datatype = 16; % single or float32
  % set slope and intercept
  nii.hdr.dime.scl_slope = 1;
  nii.hdr.dime.scl_inter = 0;

  nii.hdr.dime.dim = [4 nX nY nZ 4 1 1 1];

  nii.img = zeros(nX,nY,nZ, 3);
  nii.img(:,:,:,1) = optimalContrast;
  nii.img(:,:,:,2) = spinDensity;
  %nii.img(:,:,:,3) = -1./decayRate; % convert decay rate to time constant
  nii.img(:,:,:,3) = -decayRate; % decay rate is a better output: negative can be set to 0
  nii.img(:,:,:,4) = inverseFitError;

end



%-----------------------------------------------

function fieldOffset = internal_fitMultiEchoPhase(teList, inNiiMagnitude, inNiiPhase, carrierOffset, intensityScale, denominator, totalWeight, weightedX )
  % fit field map
  % inNiiPhase input (in radians)
  % carrierOffset input frequency offset in same time units as list of TEs
  % inNiiPhase and carrierOffset must be supplied as a pair
  % fitMask ~=0 is where to fit field map
  %
  % output is always single float
  %   field map, if phase is supplied (in supplied time units)

  assert(nargin == 8 , 'Supply filenames: teList, inMagnitudeTimeseries, inPhaseTimeseries, carrierOffset, intensityScale, denominator, totalWeight, weightedX');

  %-------------------------------
  % load magnitude data
  % We will load one volume at a time to conserve memory.
  % get matrix size from header
  %-------------------------------
  fullNiiHeader = load_nii_hdr(inNiiMagnitude); %nii.hdr

  nX  = fullNiiHeader.dime.dim(2);
  nY  = fullNiiHeader.dime.dim(3);
  nZ  = fullNiiHeader.dime.dim(4);
  nT  = fullNiiHeader.dime.dim(5);

  assert(nT >= 2, 'At least 2 time points on input are expected');

  %-------------------------------
  % load header for phase data
  %-------------------------------
  fullNiiHeader = load_nii_hdr(inNiiPhase); %nii.hdr

  nXPhase  = fullNiiHeader.dime.dim(2);
  nYPhase  = fullNiiHeader.dime.dim(3);
  nZPhase  = fullNiiHeader.dime.dim(4);
  nTPhase  = fullNiiHeader.dime.dim(5);

  assert(nXPhase == nX && nYPhase == nY && nZPhase == nZ && nTPhase == nT, 'inMagnitudeTimeseries and inPhaseTimeseries must have the same size');

  %-------------------------------
  % fit field offset
  % model: phase = slope*x + intercept
  %   slope = fieldOffset, intercept = initial phase,
  % for this to work, phase must be unwrapped. If delta TE is too long
  % and jump by more than 2pi takes place, we will not recover correct field offset.
  % To unwrap, need to compare to previous time point.
  % Wrap is when jump is larger than pi. Therefore, for the first time point, reference is zero:
  % Zero is never more than pi away from phase if phase is in (-pi, pi)
  % For others, previous unwrapped value is the reference
  %
  % S is magnitude image with noise independent of time
  % therefore, noise in phase is 1/S and weight in least squares is S^2; same as in the decay fit
  % we will use imaginary part for field to compute decay and field at the same time
  %-------------------------------

  % list of echo times (teList) can be in arbitrary order for decay fit
  % for fitting field map, we have to order by TE because phase is limited to +/-2pi
  % order in ascending order of TE
  [x, volumeOrder] = sort(teList);
  assert(length(x) == nT, 'Number of TEs differs from number of time points in NIFTI file');

  weightedY  = zeros(nX,nY,nZ); % sum of w*y
  weightedXY = zeros(nX,nY,nZ); % sum of w*x*y

  phaseReference = zeros(nX,nY,nZ);
  fieldOffset = zeros(nX,nY,nZ);

  for k = 1:nT

    % load one time point at a time
    nii = load_untouch_nii(inNiiMagnitude, volumeOrder(k));

    % Since we use load_untouch_nii, no scale/intercept is
    % applied. Apply them and convert to floating number to allow NaN
    nii.img = nii.hdr.dime.scl_slope * single(nii.img) + nii.hdr.dime.scl_inter;

    weight = (nii.img / intensityScale) .* (nii.img / intensityScale);

    %-------------------------------
    % load phase data
    %-------------------------------
    % load one time point at a time
    nii = load_untouch_nii(inNiiPhase, volumeOrder(k));
    nii.img = nii.hdr.dime.scl_slope * single(nii.img) + nii.hdr.dime.scl_inter;

    % unwrap phase by making sure new phase is within (-pi,pi) of the previous
    % offset center frequency to carrierOffset
    nii.img = nii.img - phaseReference - 2*pi*x(k) * carrierOffset;
    nii.img = phaseReference + atan2(sin(nii.img), cos(nii.img));
    phaseReference = nii.img;

    % elements of weighted least squares

    weightedY  = weightedY  + weight.*nii.img;
    weightedXY = weightedXY + weight.*x(k).*nii.img;

  end; % for k

  %-------------------------------
  % compute and output field map
  %-------------------------------
  fieldOffset = (totalWeight.*weightedXY - weightedX.*weightedY) ./ denominator; % B0 per time unit
  fieldOffset = fieldOffset/2/pi; % from radians per time unit to revolutions per time unit
  % offset center frequency to carrierOffset

  fieldOffset = fieldOffset + carrierOffset;

end


