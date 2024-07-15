function senmap_run3D(senmapdir,read_shift,phase_shift)
% SENMAP_RUN(senmapdir)
%   reconstructs the data collected with the "senmap" sequence and 
%   calculates the coil sensitivity and magnetic field maps
%   
%   INPUTS:
%   - senmapdir:    dirname of senmap data

cd(senmapdir)

recoInfo = setup_recon;
% recoInfo.nsl = nsl;
recoInfo.shft_FOVy = phase_shift;
recoInfo.shft_FOVx = read_shift;
recoInfo.num_adc = 1;
recoInfo.triple = 0;
save recoInfo.mat recoInfo

spiral_recon(recoInfo);

% create_sen_map_rough(recoInfo);
create_sen_map_mask(recoInfo);

create_sos_mask(recoInfo);

recoInfo.rec_type = 7;
spiral_recon(recoInfo);

create_field_map(recoInfo);

cd ..