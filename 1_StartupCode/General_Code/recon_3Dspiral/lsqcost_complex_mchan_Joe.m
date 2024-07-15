function thecost = lsqcost_complex_mchan_Joe(c_vec, meas_img, orig_img, mask)
    %change_log
    % December 12, 2011: JLH - The cost function was rewritten to reflect
    %   equation 3 in the Van 2011 TMI paper.  It formerly minimized
    %   equation 4.
    % December 12, 2011: JLH - the cost has real and imaginary components,
    %   instead of magnitude.

    % fitting on the real part of the signal model
    [nx ny nz nc] = size(meas_img);
    xpts = -1/2:1/nx:(1/2 - 1/nx);  % x coordinate
    ypts = -1/2:1/ny:(1/2 - 1/ny);  % y coordinate
    zpts = -1/2:1/nz:(1/2 - 1/nz);
    xmtx = repmat(col(xpts), [1, ny, nz]);
    ymtx = repmat(col(ypts).', [nx, 1, nz]);
    zmtx = repmat(col(zpts).', [nx*ny 1]);
    zmtx = reshape(zmtx, [nx ny nz]);
%     slice_weight = ones(size(meas_img));
%     for ii = 1:nz
%         slice_weight(:,:,ii,:) = abs(ii-(nz+1)/2);
%     end

    add_phase = c_vec(1)*xmtx + c_vec(2)*ymtx + c_vec(3)*zmtx +c_vec(4);
    clear xmtx ymtx zmtx xpts ypts zpts
    mymask = repmat(mask, [1 1 1 nc]);
    add_phase = repmat(add_phase, [1 1 1 nc]);
    model_img = 1.*exp(1j*(angle(orig_img) + add_phase));
    %resid = col(slice_weight).*abs(col((meas_img(mymask==1)))).*(col((meas_img(mymask==1))./(abs(meas_img(mymask==1))) - model_img(mymask==1)));
    resid = abs(col((meas_img(mymask==1)))).*(col((meas_img(mymask==1))./(abs(meas_img(mymask==1))) - model_img(mymask==1)))*1e12;
%     resid = (col((meas_img(mymask==1))./(abs(meas_img(mymask==1))) - model_img(mymask==1)));
    thecost = cat(1,real(resid),imag(resid));
%     thecost = abs(resid);
end