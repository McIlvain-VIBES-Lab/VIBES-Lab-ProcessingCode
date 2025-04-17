% Function to register a t2 image to a reference image, then apply the
% transform to a property stack. Property stack can have multiple property
% types (stored as a 4d stack)


function [propstack_reg,t2stack_reg]=Register_Stack_to_T2(t2stacki,propstacki,t2stackref,voxsize_mm)
    % COnvert to nifti
    
    t2stackinii = make_nii(flipdim(flipdim(permute(t2stacki,[2 1 3]),1),2),[voxsize_mm(2) voxsize_mm(1) voxsize_mm(3)]);
    save_nii(t2stackinii,'t2stacki.nii')
    t2stackrefnii = make_nii(flipdim(flipdim(permute(t2stackref,[2 1 3]),1),2),[voxsize_mm(2) voxsize_mm(1) voxsize_mm(3)]);
    save_nii(t2stackrefnii,'t2stackref.nii')
    
    !$FSLDIR/bin/flirt -in t2stacki.nii -ref t2stackref.nii -out t2_to_ref.nii -omat t2_to_ref.mat -dof 6

    !gunzip -f t2_to_ref.nii.gz
    tmp = load_nii('t2_to_ref.nii');
    t2stack_reg = double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3]));
 
    s=size(propstacki);
    propstackreg=zeros(s);
    for ii=1:s(4)
        propstackinii = make_nii(flipdim(flipdim(permute(propstacki(:,:,:,ii),[2 1 3]),1),2),[voxsize_mm(2) voxsize_mm(1) voxsize_mm(3)]);
        save_nii(propstackinii,'propstacki.nii');
        !$FSLDIR/bin/flirt -in propstacki.nii -ref t2stackref.nii -out propstack_Reg.nii -init t2_to_ref.mat -applyxfm;

        !gunzip -f propstack_Reg.nii.gz
        tmp = load_nii('propstack_Reg.nii');
        propstack_reg(:,:,:,ii) = double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3]));
    end
      
    showfigs=true;
    if(showfigs)
        montagestack(t2stackref);title('t2 stack_ref');pause(0.1);drawnow;pause(0.1);drawnow;
        montagestack(t2stack_reg);title('t2 stack registered');pause(0.1);drawnow;pause(0.1);drawnow;
        montagestack(t2stacki);title('t2 stack unregistered');pause(0.1);drawnow;pause(0.1);drawnow;
    end
end

