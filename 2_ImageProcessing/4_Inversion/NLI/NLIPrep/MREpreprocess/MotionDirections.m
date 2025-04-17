function MotionDirection(base)
dspfile=dir('*.v3c');

[dspfull]=load(dspfile(1).name); 
re_dsp=dspfull(:,[2 4 6]);
im_dsp=dspfull(:,[3 5 7]);

%cases={'A1','A2','A3','A4','A5','A6',...
%       'D1','D2','D3','D4','D5','D6',...
%       'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18',...
%       'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18'};
cases={'B6','C13','B10','C9'};
for ii=1:length(cases); 
    dirname=strcat('DirIndex',cases(ii),'.mat');
    load(dirname{1})
    DI=DirIndex(1:3,4:6); 
    re_dsp_new=re_dsp*DI;
    im_dsp_new=im_dsp*DI;
    dspfulltemp=[dspfull(:,1),re_dsp_new(:,1),im_dsp_new(:,1),re_dsp_new(:,2),im_dsp_new(:,2),re_dsp_new(:,3),im_dsp_new(:,3)];
    dsp_fnamenew=strcat(base,'_MotionPerm',cases(ii),'.dsp'); 
    fid=fopen(dsp_fnamenew{1},'w');
    fprintf(fid,'%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E \n',dspfulltemp');
    fclose(fid);    
end


for ii=1:length(cases); 
    
fid=fopen(['runfile-v8_porotet_MotionPerm',cases{ii},'.dat'],'w');
fprintf(fid,'Subzone Reconstruction Data File \n');
fprintf(fid,'Problem Type <0=Forward Problem Only> <1+=Inverse Problem> \n');
fprintf(fid,'1 \n');
fprintf(fid,'Node File: \n');
fprintf(fid,[base,'.nod\n']);
fprintf(fid,'Element File: \n');
fprintf(fid,[base,'.elm\n']);
fprintf(fid,'Boundary Node File: \n');
fprintf(fid,[base,'.bnod\n']);
fprintf(fid,'Pressure Boundary condition file (format: <#nbc, bnod, typ, realval, imagval>) \n');
fprintf(fid,[base,'.pbcs\n']);
fprintf(fid,'Region Stack File: \n');
fprintf(fid,[base,'.regstack\n']);
fprintf(fid,'Measured Displacement File: \n');
fprintf(fid,'1 \n');
fprintf(fid,'    1.3300d0\n');
fprintf(fid,[base,'_MotionPerm',cases{ii},'.dsp \n']);
fprintf(fid,'Initial Solution File: \n');
fprintf(fid,'0 \n');
fprintf(fid,'Output File Stem: \n');
fprintf(fid,['/ihome/lsolamen/data/data_space/TofuGelInclusionSet3/',base,'/MotionPerm',cases{ii},'/MotionPerm',cases{ii},'\n']);
fprintf(fid,'Print Detailed Runtime Debugging and Execution Comments <verb> <file>:\n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'Material Model <1=isotropic> <2=orthotropic> <3=iso compress>  <4=poro>: \n');
fprintf(fid,'4 \n');
fprintf(fid,'Number of Material Properties: \n');
fprintf(fid,'4 \n');
fprintf(fid,'Material Description Style <1=nodal> <2=element> [<shear modulus> <lambda modulus> <log10(HC)> <porosity>]: \n');
fprintf(fid,'1,1,1,1 \n');
fprintf(fid,'Number of Parameters per Material Point: \n');
fprintf(fid,'1,1,1,1 \n');
fprintf(fid,'Property Scalars: <Shear Mod, Lamba Mod, Hydraulic Conductivity, Porosity> \n');
fprintf(fid,' 3300.d0, 1188.d0,  4300.000d0,     0.000d0,  -8.00000d0, -24.00000d0,  0.200d0,  0.000d0\n');
fprintf(fid,'Other Poroealstic constants assumed constant <rho, rhof, rhoa> \n');
fprintf(fid,'1020.d0, 1000.d0, 150.d0 \n');
fprintf(fid,'Multifrequency Indicator \n');
fprintf(fid,'.false.,.false.,.false.,.false. \n');
fprintf(fid,'Multifrequency Alpha \n');
fprintf(fid,'0.d0,0.d0,0.d0,0.d0 \n');
fprintf(fid,'0.d0,0.d0,0.d0,0.d0 \n');
fprintf(fid,'Number of Material Property Mesh Resolutions: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Number of material mesh structures \n');
fprintf(fid,'3 \n');
fprintf(fid,'Material mesh iteration stucture (one less than number of structures)  \n');
fprintf(fid,'3,6 \n');
fprintf(fid,'Material Property Mesh Resolutions (x,y,z): \n');
fprintf(fid,'0.0030000 , 0.0030000 , 0.0030000 \n');
fprintf(fid,'0.0200000 , 0.0200000 , 0.0200000 \n');
fprintf(fid,'0.0400000 , 0.0400000 , 0.0400000 \n');
fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n');
fprintf(fid,' 1 1 1 3 \n');
fprintf(fid,'3 3 3 3 \n');
fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n');
fprintf(fid,' 1 1 1 3 \n');
fprintf(fid,'3 3 3 3 \n');
fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n');
fprintf(fid,' 1 1 1 3 \n');
fprintf(fid,'3 3 3 3 \n');
fprintf(fid,'Reconstruction Indicators: \n');
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'Property Estimate Variance Calculation: \n');
fprintf(fid,'.false. \n');
fprintf(fid,'Zone Sizes (not including overlap factor, actual size = (1+2*ovlp)*siz) [x y z]: \n');
fprintf(fid,' 1.95633e-02, 1.95633e-02, 1.95633e-02\n');
fprintf(fid,'Zone Overlap Percentage [x y z]: \n');
fprintf(fid,'0.150d0,0.150d0,0.150d0\n');
fprintf(fid,'Iteration Limits [global zone]: \n');
fprintf(fid,'750,1 \n');
fprintf(fid,'Minimum Parameter Update Size [global zone line]: \n');
fprintf(fid,'0.1d0, 0.1d0, 0.1d0\n');
fprintf(fid,'Number of Zone Iteration Structures: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Iteration Limits for Zone Iteration Structures (NOTE:  provide one less limit than # of structures!!!): \n');
fprintf(fid,'10,150 \n');
fprintf(fid,'Zone Iteration Structures [<# of CG iters> <# of GN iters> <!!! QN CURRENTLY UNAVAILABLE !!!> <# of line search iters>]: \n');
fprintf(fid,'0,0,1,1 \n');
fprintf(fid,'0,0,2,2 \n');
fprintf(fid,'0,0,3,2 \n');
fprintf(fid,'Number of Processors per MUMPS Communicator \n');
fprintf(fid,'1 \n');
fprintf(fid,'Maximum Amount of RAM per Processor [MB] \n');
fprintf(fid,'999999999 \n');
fprintf(fid,'ooooooooooo  REGULARIZATION DESCRIPTORS  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n');
fprintf(fid,'Regularization Indicators [<TV> <SF> <TK> <MQ> <JW> <constraints> <CG residual> <VH> <McG> <soft prior>]: \n');
fprintf(fid,'.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true.,.true.,.false. \n');
fprintf(fid,'Number of constant regularization iterations (final N iterations use final regularization weights) \n');
fprintf(fid,'30 \n');
fprintf(fid,'Number of Parameters to Treat with Total Variation: \n');
fprintf(fid,'4 \n');
fprintf(fid,'TV Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'TV Delta Values: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-19,1.d-19,1.d-19,1.d-19 \n');
fprintf(fid,'1.d-19,1.d-19,1.d-19,1.d-19 \n');
fprintf(fid,'TV Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'TV Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'Number of Parameters to Treat with Spatial Filtering: \n');
fprintf(fid,'4 \n');
fprintf(fid,'SF Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'SF Gaussian filter initial Widths: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'SF Final Gaussian width: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'Number of Parameters to Treat with Tikhonov Regularization: \n');
fprintf(fid,'4 \n');
fprintf(fid,'TK Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'TK Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'TK Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'Number of Parameters to Treat with Marquardt Regularization: \n');
fprintf(fid,'4 \n');
fprintf(fid,'Distance of alpha from 1 at which MQ weights are adjusted: \n');
fprintf(fid,'0.25d0 \n');
fprintf(fid,'MQ Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'MQ Weight Delta: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'MQ Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d2,1.d2,1.d2,1.d2 \n');
fprintf(fid,'1.d2,1.d2,1.d2,1.d2 \n');
fprintf(fid,'MQ Minimum Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'Number of Parameters to Treat with Joachimowicz Regularization: \n');
fprintf(fid,'4 \n');
fprintf(fid,'JW Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'JW Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0 \n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0 \n');
fprintf(fid,'JW Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0 \n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0 \n');
fprintf(fid,'Number of Parameters to Treat with Constraints: \n');
fprintf(fid,'4 \n');
fprintf(fid,'Constraint Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'Constraint Weights <converges to exact solution as weight --> inf>: \n');
fprintf(fid,'1.d-14,1.d-14,1.d-14,1.d-14 \n');
fprintf(fid,'Constraint Minimums: \n');
fprintf(fid,'200.d0, 0.d0, 200.d0, 0.d0, -11.d0, -24.d0, 0.01d0, 0.d0  \n');
fprintf(fid,'Constraint Maximums: \n');
fprintf(fid,'2.d4, 5.d5, 5.d4, 5.d5, -5.d0, -4.d0, 0.8d0, 0.8d0 \n');
fprintf(fid,'CG Residual Scaling Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'CG Residual Scaling Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'Van Houten Regularization Level: \n');
fprintf(fid,'1.2d0 \n');
fprintf(fid,'Number of parameters to treat with soft prior: \n');
fprintf(fid,'4 \n');
fprintf(fid,'SP Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'SP Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-17,1.d-17,1.d-17,1.d-17 \n');
fprintf(fid,'1.d-17,1.d-17,1.d-17,1.d-17 \n');
fprintf(fid,'SP Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-17,1.d-17,1.d-17,1.d-17 \n');
fprintf(fid,'1.d-17,1.d-17,1.d-17,1.d-17 \n');
fprintf(fid,'SP start iteration delay (first N iterations do not use soft prior \n');
fprintf(fid,'0 \n');
fclose(fid);

fid2=fopen(['MPICH2-Submitv8_poro',cases{ii}], 'w');
fprintf(fid2, '#!/bin/bash -l \n');
fprintf(fid2, '# declare a name for this job \n');
fprintf(fid2, ['#PBS -N',base,'_',cases{ii},'\n']);
fprintf(fid2, '# request the queue (enter the possible names, if omitted, serial is the default) \n');
fprintf(fid2, '#PBS -q default  \n');
fprintf(fid2, '# request  node  \n');
fprintf(fid2, '#PBS -l nodes=4:ppn=8 \n');
fprintf(fid2, '#PBS -l feature=intel \n');
fprintf(fid2, '# request some hours of wall time  \n');
fprintf(fid2, '#PBS -l walltime=36:00:00  \n');
fprintf(fid2, '#combine PBS standard output and error files  \n');
fprintf(fid2, '#PBS -j oe  \n');
fprintf(fid2, '# mail is sent to you when the job starts and when it terminates or aborts  \n');
fprintf(fid2, '##PBS -m bea  \n');
fprintf(fid2, '# specify your email address  \n');
fprintf(fid2, '##PBS -M matthew.d.mcgarry@dartmouth.edu  \n');
fprintf(fid2, '# Change to Submission Directory  \n');
fprintf(fid2, 'cd $PBS_O_WORKDIR  \n');
fprintf(fid2, '# run the program \n');
fprintf(fid2, 'cat $PBS_NODEFILE | uniq > node_file \n');
fprintf(fid2, 'nnodes=$(cat node_file | wc -l) \n');
fprintf(fid2, 'nprocs=$(cat $PBS_NODEFILE | wc -l) \n');
fprintf(fid2, 'export MKL_NUM_THREADS=1 \n');
fprintf(fid2, 'echo Nodes $nnodes \n');
fprintf(fid2, 'echo Procs $nprocs \n');
fprintf(fid2, ['mpirun -np $nprocs /ihome/lsolamen/code/MREv8/MRE-Zone.mp2 runfile-v8_porotet_MotionPerm',cases{ii},'.dat \n']);
fprintf(fid2, 'exit 0  \n');
fclose(fid2);

clear fid fid2
end


end

   
% ADirIndexnames=dir('DirIndexA*');
% BDirIndexnames=dir('DirIndexB*');
% CDirIndexnames=dir('DirIndexC*');
% DDirIndexnames=dir('DirIndexD*');
% clear DirIndex;
% 
% 
% 
% for ii=1:length(ADirIndexnames);
%     load(ADirIndexnames(ii).name); 
%     DI=DirIndex(1:3,4:6); 
%     re_dsp_new=re_dsp*DI;
%     im_dsp_new=im_dsp*DI;
%     dspfulltemp=[dspfull(:,1),re_dsp_new(:,1),im_dsp_new(:,1),re_dsp_new(:,2),im_dsp_new(:,2),re_dsp_new(:,3),im_dsp_new(:,3)];
%     dsp_fnamenew=strcat('set2_MotionPermA',num2str(ii),'.dsp'); 
%     fid=fopen(dsp_fnamenew,'w');
%     fprintf(fid,'%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E \n',dspfulltemp');
%     fclose(fid);
% end
% 
% clear DirIndex;
% 
% 
% for ii=1:length(BDirIndexnames);
%     load(BDirIndexnames(ii).name); 
%     DI=DirIndex(1:3,4:6); 
%     re_dsp_new=re_dsp*DI;
%     im_dsp_new=im_dsp*DI;
%     dspfulltemp=[dspfull(:,1),re_dsp_new(:,1),im_dsp_new(:,1),re_dsp_new(:,2),im_dsp_new(:,2),re_dsp_new(:,3),im_dsp_new(:,3)];
%     dsp_fnamenew=strcat('set2_MotionPermB',num2str(ii),'.dsp');  
%     fid=fopen(dsp_fnamenew,'w');
%     fprintf(fid,'%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E \n',dspfulltemp');
%     fclose(fid);
% end
% 
% clear DirIndex;
% 
% 
% for ii=1:length(CDirIndexnames);
%     load(CDirIndexnames(ii).name); 
%     DI=DirIndex(1:3,4:6); 
%     re_dsp_new=re_dsp*DI;
%     im_dsp_new=im_dsp*DI;
%     dspfulltemp=[dspfull(:,1),re_dsp_new(:,1),im_dsp_new(:,1),re_dsp_new(:,2),im_dsp_new(:,2),re_dsp_new(:,3),im_dsp_new(:,3)];
%     dsp_fnamenew=strcat('set2_MotionPermC',num2str(ii),'.dsp'); 
%     fid=fopen(dsp_fnamenew,'w');
%     fprintf(fid,'%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E \n',dspfulltemp');
%     fclose(fid);
% end
% 
% clear DirIndex;
% 
% 
% for ii=1:length(DDirIndexnames);
%     load(DDirIndexnames(ii).name); 
%     DI=DirIndex(1:3,4:6); 
%     re_dsp_new=re_dsp*DI;
%     im_dsp_new=im_dsp*DI;
%     dspfulltemp=[dspfull(:,1),re_dsp_new(:,1),im_dsp_new(:,1),re_dsp_new(:,2),im_dsp_new(:,2),re_dsp_new(:,3),im_dsp_new(:,3)];
%     dsp_fnamenew=strcat('set2_MotionPermD',num2str(ii),'.dsp'); 
%     fid=fopen(dsp_fnamenew,'w');
%     fprintf(fid,'%6d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E \n',dspfulltemp');
%     fclose(fid);
% end


   %cases={'B15'}

