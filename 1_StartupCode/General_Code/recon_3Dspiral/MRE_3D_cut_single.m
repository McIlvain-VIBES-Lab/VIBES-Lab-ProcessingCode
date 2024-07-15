function MRE_3D_cut_single(cimg,kk)

% if ~exist('nImages','var')
%    nImages =25;
% elseif isempty(nImages)
%    nImages =25;
% end

%     nImages =25;
%     nImages =16;
%     nSlabs = 4;
    nSlabs=10;
%     nSlices =60;
    nSlices =8;
%     dist_fact = .20;
    dist_fact =0.25;
    szx = 2; szy = 2; szz = 2;
%     szx = 1; szy = 1; szz = 1;
%     szx = 1.25; szy = 1.25; szz = 1.25;
%     szx = 0.8; szy = 0.8; szz = 0.8;
    start_idx = nSlices*dist_fact/2;
    stop_idx = nSlices - start_idx;
%     start_idx = 0;
%     stop_idx = 8;
    slices_keep = stop_idx-start_idx;

%     for kk = 1:nImages
%     for kk = 2
%         eval(sprintf('load img_%d', kk));
%         eval(sprintf('load img_grid_%d', kk));
%         imgm = permute(cimg(:,:,:),[2,1,3]);
%         imgm = permute(img_grid(:,:,:),[2,1,3]);
%         eval(sprintf('save rKICT_vol_%02d imgm', (kk)));
%         write_ana(sprintf('aKICT_%02d', (kk)), imgm, szx, szy, szz, 'FLOAT')
%         eval(sprintf('load img_grid_%d', kk));
%         img_tmp = zeros(size(img_grid,1),size(img_grid,2),nSlabs*nSlices*dist_fact);

        for ii = 1: nSlabs
            img(:,:,((ii-1)*slices_keep+1):(slices_keep*ii)) = cimg(:,:,((ii-1)*nSlices+start_idx+1):((ii-1)*nSlices+stop_idx));
        end
        eval(sprintf('save img_cut_%i img',kk))

%         img_tmp =imgm;
%         imgm = permute(img_tmp(:,:,:),[2,1,3]);
%         write_ana(sprintf('DWI_%02d', (kk)), abs(img_tmp), szx, szy, szz, 'FLOAT')
%         write_ana(sprintf('DWI_%02d', (kk)), abs(img_tmp(:,:,1:24))*1e12, szx, szy, szz, 'FLOAT')
%     end
    
    
%     !fslmerge -t merged_data DWI_*.hdr
    
    %insert BET step
%     write_ana('mask', ones(size(img_tmp)), szx, szy, szz, 'FLOAT')
%    write_ana('mask', ones(240, 240,96), szx, szy, szz, 'FLOAT')
%     write_ana('mask', ones(240, 240,96), szx, szy, szz, 'SHORT')
%    !bet DWI_01.hdr brain -f 0.3 -m
%     !bet DWI_01.hdr brain -f 0.15  -m
%     !dtifit --data=merged_data.nii.gz --out=dti --mask=brain_mask.nii.gz --bvecs=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_1b0_15D_b700/bvecs.txt --bvals=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_1b0_15D_b700/bvals.txt
%     !dtifit --data=merged_data.nii.gz --out=dti --mask=brain_mask.nii.gz --bvecs=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_1b0_6D_b700/bvec6.txt --bvals=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections//bval6_1600.txt
%     !dtifit --data=merged_data.nii.gz --out=dti --mask=brain_mask.nii.gz --bvecs=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_2b0_15D_b700/bvecs.txt --bvals=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_2b0_15D_b700/bvals.txt
%     !dtifit --data=merged_data.nii.gz --out=dti --mask=brain_mask.nii.gz --bvecs=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_2b0_15D_b1000/bvecs.txt --bvals=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_2b0_15D_b1000/bvals.txt
%     !dtifit --data=merged_data.nii.gz --out=dti --mask=brain_mask.nii.gz --bvecs=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_1b0_6D_b700/bvec6.txt --bvals=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/3D_1b0_6D_b700/bval6.txt
%  !dtifit --data=merged_data.nii.gz --out=dti --mask=brain_mask.nii.gz --bvecs=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/FSL_30D_2B0_axial/bvecs_30.txt --bvals=/home/UIUC/holtrop1/MRFIL/archive/DTIdirections/FSL_30D_2B0_axial/bvals_30.txt
end