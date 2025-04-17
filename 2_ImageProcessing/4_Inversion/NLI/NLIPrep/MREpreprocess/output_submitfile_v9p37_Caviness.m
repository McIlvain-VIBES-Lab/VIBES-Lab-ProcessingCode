function output_submitfile_v9p37_Caviness(rfilestem,rfile,ifile,outstm,inpath,outpath,freqHz,vox,modelsiz,globitr)
%output runfiles in new general format to hopefully make things simpler
%
% 
% 1 isoincomp
% 2 isoincomp_SPon
% 3 isocomp
% 4 poro_noHC
% 5 poro_noHC_SPon
% 6 poro_reconHC
% 7 poro_reconHC_SPon
% 8 tetvisc
% 9 tetvisc_SPon
% 10 NITI_anisodamp_SPoff
% 11 NITI_isodamp_SPoff
% 12 NITI_anisodamp_SPoff
% 13 NITI_isodamp_SPoff
% Parmeters to be changed

% Example Caviness Submit File Format
% #!/bin/bash -l 
% # declare a name for this job 
% #SBATCH --job-name=NLI_200807_GM_Rank8_Ep2d_dir_space_voxelmesh_G3300.v7.4.inv.iso.incomp.visc_SPoff_SF0p0015_CG2p2
% # request  nodes  
% #SBATCH --nodes=4
% #SBATCH --tasks-per-node=18
% # request memory per node  
% #SBATCH --mem=60G
% # one task per CPU  
% #SBATCH --cpus-per-task=1
% # request some hours of wall time <day-hour:min:sec> 
% #SBATCH --time=2-12:00:00  
% # do some other stuff  
% #SBATCH --partition=standard
% #SBATCH --export=NONE
% vpkg_require mumps
% . /opt/shared/slurm/templates/libexec/openmpi.sh
% ${UD_MPIRUN} /work/clj/MREcode/MREv7p40/bin/mre-zone runfile-v7p3visc.dat
% exit $mpi_rc 

procpernode=18;
if((ifile==1)||(ifile==2)) % isotropic visco
    model=1;
    Nnode=4; % Number of discovery nodes to use
    zonetim2mm=0.01; % Time in hours for 1 zone @2mm resolution to run
elseif(ifile==4)||(ifile==5)||(ifile==6)||(ifile==7) % poroelastic
    model=4;
    Nnode=4;
    zonetim2mm=0.02; 
elseif((ifile==3)||(ifile==8)||(ifile==9)) % isotopic compressible viscoelasticv
    model=3;
    Nnode=4;
    zonetim2mm=0.015; 
elseif((ifile==10)||(ifile==11)) % aniso: NITI aniso damping 
    model=5;
    Nnode=8;
    zonetim2mm=0.05; 
elseif(ifile==12)||(ifile==13) % aniso: NITI iso damping
    model=6;
    Nnode=8;
    zonetim2mm=0.05; 
else
    disp('WARNING: Nprop not defined in submitfile');
end
    
%% Estimate runtime to set wall time limit
resfac=(0.002/mean(vox))^4.6; % Runtime ^1.6 for Number of fwd prob nodes, ^3 for nodes per zone
zonesiz=0.02; % APproz zone size
nzone=ceil(modelsiz(1)/zonesiz)*ceil(modelsiz(2)/zonesiz)*ceil(modelsiz(3)/zonesiz);
runtime=ceil(1.5*zonetim2mm*resfac*globitr*ceil(nzone/(Nnode*procpernode-1))); % hours

Ndisp=length(freqHz);

% #!/bin/bash -l 
% # declare a name for this job 
% #SBATCH --job-name=NLI_200807_GM_Rank8_Ep2d_dir_space_voxelmesh_G3300.v7.4.inv.iso.incomp.visc_SPoff_SF0p0015_CG2p2
% # request  nodes  
% #SBATCH --nodes=4
% #SBATCH --tasks-per-node=18
% # request memory per node  
% #SBATCH --mem=60G
% # one task per CPU  
% #SBATCH --cpus-per-task=1
% # request some hours of wall time <day-hour:min:sec> 
% #SBATCH --time=2-12:00:00  
% # do some other stuff  
% #SBATCH --partition=standard
% #SBATCH --export=NONE
% vpkg_require mumps
% . /opt/shared/slurm/templates/libexec/openmpi.sh
% ${UD_MPIRUN} /work/clj/MREcode/MREv7p40/bin/mre-zone runfile-v7p3visc.dat
% exit $mpi_rc 

% Build MRE-submit file
        
        outstm=[outstm '.MREv9p37' rfile.ID];
        fid=fopen([inpath 'CAVSUB_' rfilestem rfile.ID],'w');
        fprintf(fid,'#!/bin/bash -l \n');
        fprintf(fid,'# declare a name for this job \n'); 
        fprintf(fid,'#SBATCH --job-name=j');fprintf(fid,'%c',outstm);fprintf(fid,'\n');
        fprintf(fid,'# request some nodes  \n');
        fprintf(fid,'#SBATCH --nodes=%i \n',Nnode);
        fprintf(fid,'#SBATCH --tasks-per-node==%i \n',procpernode);
        fprintf(fid,'# request memory per node   \n');
        fprintf(fid,'#SBATCH --mem=60G \n');        
        fprintf(fid,'# one task per CPU  \n');
        fprintf(fid,'#SBATCH --cpus-per-task=1 \n');        
        fprintf(fid,'# request some hours of wall time  \n');
        fprintf(fid,'#SBATCH --time=%i:00:00  \n',runtime);
        fprintf(fid,'# do some other stuff  \n');
        fprintf(fid,'#SBATCH --partition=standard \n');
        fprintf(fid,'#SBATCH --export=NONE \n');
        fprintf(fid,'vpkg_require mumps \n');        
        fprintf(fid,'. /opt/shared/slurm/templates/libexec/openmpi.sh \n');
        fprintf(fid,'# run the program \n');        
        fprintf(fid,'rfile=');fprintf(fid,'%c',[rfilestem rfile.ID]);fprintf(fid,'\n');
        fprintf(fid,'${UD_MPIRUN} /work/clj/MREcode/MREv9p37/bin/mre-zone $rfile \n');
        fprintf(fid,'## auto post processing \n'); % this might take a bit of messing around to get to work. 
        fprintf(fid,'srchtag=*RE..*0001*prop.01.mtr\n');
        fprintf(fid,'outstm=$(sed ''%i!d'' $rfile)\n',21+2*(Ndisp-1)); % Get output line from runfile, each disp file adds 2 lines. 
        fprintf(fid,'/update/this/path/to/this/complied/fortran/code/MREpostprocessv9.x $outstm$srchtag');fprintf(fid,'\n');
        fprintf(fid,'exit $mpi_rc   \n');
        fclose(fid);
    end
