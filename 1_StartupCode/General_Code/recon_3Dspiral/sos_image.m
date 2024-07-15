function imgsos = sos_image(img_in);
%function imgsos = sos_image(img_in);

    
    nd = ndims(img_in);
    dims = size(img_in);
    
    
    ny = dims(1);
    nx = dims(2);
    if (nd == 3)
        nz = 1;
        nc = dims(3);
    elseif (nd == 4)
        nz = dims(3);
        nc = dims(4);
    else
        sprintf('sos_image does not work for arrays of higher than 3D images')
        keyboard
    end
    
    img_in = reshape(img_in,ny,nx,nz,nc);
    
    imgsos = zeros(nx,ny,nz);
    for ii = 1:nc
        imgsos = imgsos + (abs(img_in(:,:,:,ii)).^2);
    end
    
    imgsos = sqrt(imgsos);
    
    
    %get back to input size
    imgsos = squeeze(imgsos);
