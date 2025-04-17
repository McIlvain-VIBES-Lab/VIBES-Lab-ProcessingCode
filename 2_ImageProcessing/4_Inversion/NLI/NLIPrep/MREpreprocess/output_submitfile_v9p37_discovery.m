function output_submitfile_v9p37_discovery(rfilestem,rfile,ifile,outstm,inpath,outpath,freqHz,vox,modelsiz,globitr)
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


procpernode=8;
if((ifile==1)||(ifile==2))
    model=1;
    Nnode=4; % Number of discovery nodes to use
    zonetim2mm=0.01; % Time in hours for 1 zone @2mm resolution to run
elseif(ifile==4)||(ifile==5)||(ifile==6)||(ifile==7)
    model=4;
    Nnode=4;
    zonetim2mm=0.02; 
elseif((ifile==3)||(ifile==8)||(ifile==9))
    model=3;
    Nnode=4;
    zonetim2mm=0.015; 
elseif((ifile==10)||(ifile==11))
    model=5;
    Nnode=8;
    zonetim2mm=0.05; 
elseif(ifile==12)||(ifile==13)
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
runtime=ceil(1.5*zonetim2mm*resfac*globitr*ceil(nzone/(Nnode*procpernode-1)));

Ndisp=length(freqHz);

% Build MRE-submit file
        
        outstm=[outstm '.MREv9p37' rfile.ID];
        fid=fopen([inpath 'DSUB_' rfilestem rfile.ID],'w');
        fprintf(fid,'#!/bin/bash -l \n');
        fprintf(fid,'# declare a name for this job \n'); 
        fprintf(fid,'#PBS -N ');fprintf(fid,'%c',outstm);fprintf(fid,'\n');
        fprintf(fid,'# request the queue (enter the possible names, if omitted, serial is the default) \n');
        fprintf(fid,'#PBS -q default  \n');
        fprintf(fid,'# request  node  \n');
        fprintf(fid,'#PBS -l nodes=%i:ppn=%i \n',Nnode,procpernode);
        fprintf(fid,'##PBS -l feature=cellm|cellk|cellj|cellh \n');
        fprintf(fid,'# request some hours of wall time  \n');
        fprintf(fid,'#PBS -l walltime=%i:00:00  \n',runtime);
        fprintf(fid,'#combine PBS standard output and error files  \n');
        fprintf(fid,'#PBS -j oe  \n');
        fprintf(fid,'# mail is sent to you when the job starts and when it terminates or aborts  \n');
        fprintf(fid,'##PBS -m bea  \n');
        fprintf(fid,'# specify your email address  \n');
        fprintf(fid,'##PBS -M matthew.d.mcgarry@dartmouth.edu  \n');
        fprintf(fid,'# Change to Submission Directory  \n');
        fprintf(fid,'cd $PBS_O_WORKDIR  \n');
        fprintf(fid,'# run the program \n');
        fprintf(fid,'cat $PBS_NODEFILE | uniq > node_file \n');
        fprintf(fid,'nnodes=$(cat node_file | wc -l) \n');
        fprintf(fid,'nprocs=$(cat $PBS_NODEFILE | wc -l) \n');
        fprintf(fid,'export MKL_NUM_THREADS=1 \n');
        fprintf(fid,'export OMP_NUM_THREADS=1 \n');
        fprintf(fid,'echo Nodes $nnodes \n');
        fprintf(fid,'echo Procs $nprocs \n');
		fprintf(fid,'cat $PBS_NODEFILE | uniq \n');
        %fprintf(fid,'mpirun -np $nprocs -hostfile $PBS_NODEFILE /ihome/mmcgarry/code/MRE_mumpsv5/MREv9/MRE-Zone.discov ');fprintf(fid,'%c',viscrfilev9_SP);fprintf(fid,'\n');
        fprintf(fid,'rfile=');fprintf(fid,'%c',[rfilestem rfile.ID]);fprintf(fid,'\n');
        fprintf(fid,'mpirun -np $nprocs -hostfile $PBS_NODEFILE /dartfs-hpc/rc/home/l/d29143l/Discovery_home/code/MRE_mumpsv5/MREv9p37MF/MRE-Zone.discov $rfile \n');
        fprintf(fid,'## auto post processing \n');
        fprintf(fid,'srchtag=*RE..*0001*prop.01.mtr\n');
        fprintf(fid,'outstm=$(sed ''%i!d'' $rfile)\n',21+2*(Ndisp-1)); % Get output line from runfile, each disp file adds 2 lines. 
        fprintf(fid,'/dartfs-hpc/rc/home/l/d29143l/Discovery_home/code/convcode_old/MREpostprocessv8.x $outstm$srchtag');fprintf(fid,'\n');
        fprintf(fid,'exit 0  \n');
        fclose(fid);
    end
