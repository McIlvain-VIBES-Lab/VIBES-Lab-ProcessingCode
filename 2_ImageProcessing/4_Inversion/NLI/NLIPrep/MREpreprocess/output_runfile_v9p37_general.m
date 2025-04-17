function [globitr]=output_runfile_v9p37_general(rfilestem,rfile,ifile,outstm,vox,muest,rhoest,DR,inpath,nodhmgf,elmhmgf,bcoutf,pbcf,regoutf,freqHz,dspoutf,outpath,znedgelength,znovlp)
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
% 12 NITI_anisodamp_SPon
% 13 NITI_isodamp_SPon
% Parmeters to be changed
vincomp=0.499; % Default Incompressible Poisson Ratio
vcomp=0.4; % Default compressible Poisson ratio
vporo=0.25; % v=0.25 means G=L


if((ifile==1)||(ifile==2)||(ifile==3)||(ifile==8)||(ifile==9))
    nprop=3;
    if(ifile<8)
        globitr=100;
    else
        globitr=300;
    end        
elseif((ifile==4)||(ifile==5)||(ifile==6)||(ifile==7))
    nprop=4;
    globitr=300;
elseif((ifile==10)||(ifile==11)||(ifile==12)||(ifile==13))
    nprop=5;
    globitr=100;
else
    disp('WARNING: Nprop not defined in runfile');
end
    
    

Ndisp=length(freqHz);

        outstm=[outstm  '_G' sprintf('%0.0f',muest) '.v9p37' rfile.ID];
        fid=fopen([inpath rfilestem rfile.ID],'w');
        fprintf(fid,'Subzone Reconstruction Data File \n');
        fprintf(fid,'Problem Type <0=Forward Problem Only> <1+=Inverse Problem> \n');
        fprintf(fid,'1 \n');
        fprintf(fid,'Node File: \n');
        fprintf(fid,'%c',nodhmgf);fprintf(fid,'\n');
        fprintf(fid,'Element File: \n');
        fprintf(fid,'%c',elmhmgf);fprintf(fid,'\n');
        fprintf(fid,'Boundary Node File: \n');
        fprintf(fid,'%c',bcoutf);fprintf(fid,'\n');
        fprintf(fid,'Pressure BC file: \n');
        if((ifile==4)||(ifile==5)||(ifile==6)||(ifile==7))
            fprintf(fid,'%c',pbcf);fprintf(fid,'\n');
        else
            fprintf(fid,'UNUSED');fprintf(fid,'\n');
        end
        fprintf(fid,'Region Stack File: \n');
        fprintf(fid,'%c',regoutf);fprintf(fid,'\n');
        fprintf(fid,'Measured Displacement File: \n');
        fprintf(fid,'%i',Ndisp);fprintf(fid,'\n');
        for idisp=1:Ndisp
            fprintf(fid,'%10.4fd0',freqHz(idisp));fprintf(fid,'\n');
            fprintf(fid,'%c',dspoutf{idisp});fprintf(fid,'\n');
        end
        fprintf(fid,'Initial Solution File: \n');
        fprintf(fid,'0 \n');
        fprintf(fid,'Output File Stem: \n');
        fprintf(fid,'%c',outpath);fprintf(fid,'%c',outstm);fprintf(fid,'\n');
        fprintf(fid,'Print Detailed Runtime Debugging and Execution Comments <verb> <file>:\n');
        fprintf(fid,'.false.,.false. \n');
        fprintf(fid,'Material Model <1=isotropic> <2=orthotropic> <3=iso compress> <4=poro> <5=NITIanisodamp> <6=NITIisodamp: \n');
        if((ifile==1)||(ifile==2))
            fprintf(fid,'1 \n');
        elseif((ifile==3)||(ifile==8)||(ifile==9))
            fprintf(fid,'3 \n');
        elseif((ifile==4)||(ifile==5)||(ifile==6)||(ifile==7))
            fprintf(fid,'4 \n');
        elseif((ifile==10)||(ifile==12))
            fprintf(fid,'5 \n');
        elseif((ifile==11)||(ifile==13))
            fprintf(fid,'6 \n') ;      
        else
            disp('WARNING: MODEL NOT DEFINED IN RUNFILE')
        end
        fprintf(fid,'Number of Materia)l Properties: \n');
        fprintf(fid,[int2str(nprop) '\n']);
        fprintf(fid,'Material Description Style <1=nodal> <2=element> [<shear modulus> <density> <bulk modulus>]: \n');
        str=sprintf('%c',repmat('1,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Number of Parameters per Material Point: \n');
        str=sprintf('%c',repmat('1,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Property Scalars: \n');
        if((ifile==1)||(ifile==2))
            fprintf(fid,'%5.0f.d0,%5.0f.d0,%5.0f.d0,%8.5fd0,%10.0f.d0,%9.5fd0',[muest 2*DR*muest rhoest -0.001 2*muest*(1+vincomp)/(3*(1-2*vincomp)) 0.0001]);fprintf(fid,'\n');
        elseif((ifile==3)||(ifile==8)||(ifile==9))
                fprintf(fid,'%5.0f.d0,%5.0f.d0,%5.0f.d0,%8.5fd0,%10.0f.d0,%11.7fd0',[muest 2*DR*muest rhoest -0.001 2*muest*vcomp/(1-2*vcomp) 0.0000001]);fprintf(fid,'\n');
        elseif((ifile==4)||(ifile==5)||(ifile==6)||(ifile==7))
                fprintf(fid,'%5.0f.d0,%11.7fd0,%5.0f.d0,%11.7fd0,%8.5fd0,%9.5fd0,%9.5fd0,%11.7fd0',[muest 0.0001 2*muest*vporo/(1-2*vporo) 0.0001 -8 -50 0.2 0.000001]);fprintf(fid,'\n');
                fprintf(fid,'Other Poroealstic constants assumed constant <rho, rhof, rhoa> \n');
                fprintf(fid,'1020.d0, 1000.d0, 150.d0 \n');
        elseif((ifile==10)||(ifile==12))
                fprintf(fid,'%5.0f.d0,%5.0f.d0,%5.0f.d0,%5.0f.d0,%5.0f.d0,%5.0f.d0,%.1fd0,%.5fd0,%9.5fd0,%.5fd0',[muest 2*DR*muest muest 2*DR*muest 2*muest*(1+vincomp) 2*DR*2*muest*(1+vincomp) rhoest -0.0001 2*muest*(1+vincomp)/(3*(1-2*vincomp)) 0.0001]);fprintf(fid,'\n');
        elseif((ifile==11)||(ifile==13))
                fprintf(fid,'%5.0f.d0,%5.0f.d0,%.3fd0,%.3fd0,%.3fd0,%.3fd0,%.1fd0,%.5fd0,%9.5fd0,%.5fd0',[muest 2*DR*muest 0 0 0 0 rhoest -0.0001 2*muest*(1+vincomp)/(3*(1-2*vincomp)) 0.0001]);fprintf(fid,'\n');
        end
        if(ifile==10)||(ifile==11)||(ifile==12)||(ifile==13) % DTI info
            fprintf(fid,'Fiber Direction Specification <0 = x direction, 1=by region, 2=regstack file> \n');
            fprintf(fid,'2 \n');
            fprintf(fid,'Fiber direction info <either 3 regstack file names or numdirs;x y z dir1;x y z dir2....> \n');
            fprintf(fid,[nodhmgf(1:end-7) 'DTIstackx \n']);
            fprintf(fid,[nodhmgf(1:end-7) 'DTIstacky \n']);
            fprintf(fid,[nodhmgf(1:end-7) 'DTIstackz \n']);
        end
        fprintf(fid,'Zone displacement boundary condition option: 1=measured data, 2=seperate file (needs disp files on the following lines),3=calculated displacements \n');
        fprintf(fid,'1 \n');
        if(ifile==4)||(ifile==5)% Poro no recon HC
            fprintf(fid,'Poroelastic zone pressure option <1:p=0, 2=pressure fwd prob, 3=global forward solve, 4=from file (needs pressure files on following lines)> \n');
            fprintf(fid,'2 \n'); 
            fprintf(fid,'Poroelastic fluid source file <.true.,.false.> Followed by filenames if true \n');
            fprintf(fid,'.false. \n');
        elseif(ifile==6)||(ifile==7) % Poro no recon HC, somethings up with the PFP and HC reconstruction. Cant figure out what.
            fprintf(fid,'Poroelastic zone pressure option <1:p=0, 2=pressure fwd prob, 3=global forward solve, 4=from file (needs pressure files on following lines)> \n');
            fprintf(fid,'1 \n'); 
            fprintf(fid,'Poroelastic fluid source file <.true.,.false.> Followed by filenames if true \n');
            fprintf(fid,'.false. \n');
        end
		fprintf(fid,'Multifrequency Indicator \n');
        str=sprintf('%c',repmat('.false.,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Multifrequency Alpha \n');
        str=sprintf('%c',repmat('0.d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('0.d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        if(ifile<10) % No fancy tricks with property mesh resolution
            fprintf(fid,'Number of Material Property Mesh Resolutions: \n');
            fprintf(fid,'2 \n');
            fprintf(fid,'Number of material mesh structures \n');
            fprintf(fid,'2 \n');
            fprintf(fid,'Material mesh iteration stucture (one less than number of structures)  \n');
            fprintf(fid,[int2str(globitr+1) '\n']);
            fprintf(fid,'Material Property Mesh Resolutions (x,y,z): \n');
            fprintf(fid,'%9.7f , %9.7f , %9.7f \n',vox);
            fprintf(fid,'%9.7f , %9.7f , %9.7f \n',10.*vox);
            if(ifile==1)||(ifile==2) % Iso incomp
                mshindr='1 2 1';
                mshindi='1 2 1';
            elseif(ifile==3)||(ifile==8)||(ifile==9) % Iso comp
                mshindr='1 2 1';
                mshindi='1 2 2';
            elseif(ifile==4)||(ifile==5) % Poro no hc recon
                mshindr='1 1 2 2';
                mshindi='2 2 2 2';
            elseif(ifile==6)||(ifile==7) % Poro hc recon
                mshindr='1 1 1 2';
                mshindi='2 2 2 2';            
            end            
            fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
            fprintf(fid,[mshindr ' \n']);
            fprintf(fid,[mshindi ' \n']);
            fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
            fprintf(fid,[mshindr ' \n']);
            fprintf(fid,[mshindi ' \n']); 
            if(ifile==1)||(ifile==2) % Iso incomp
                fprintf(fid,'Reconstruction Indicators: \n');
                fprintf(fid,'.true.,.true. \n');
                fprintf(fid,'.false.,.false. \n');
                fprintf(fid,'.false.,.false. \n');
            elseif(ifile==3)||(ifile==8)||(ifile==9) % Iso comp
                fprintf(fid,'Reconstruction Indicators: \n');
                fprintf(fid,'.true.,.true. \n');
                fprintf(fid,'.false.,.false. \n');
                fprintf(fid,'.true.,.false. \n');
            elseif(ifile==4)||(ifile==5) % Poro no hc recon
                fprintf(fid,'Reconstruction Indicators: \n');
                fprintf(fid,'.true.,.false. \n');
                fprintf(fid,'.true.,.false. \n');
                fprintf(fid,'.false.,.false. \n');
                fprintf(fid,'.false.,.false. \n');
            elseif(ifile==6)||(ifile==7) % Poro hc recon
                fprintf(fid,'Reconstruction Indicators: \n');
                fprintf(fid,'.true.,.false. \n');
                fprintf(fid,'.false.,.false. \n');
                fprintf(fid,'.true.,.false. \n'); 
                fprintf(fid,'.false.,.false. \n');
            end
        else % Change mesh resolution for aniso inversions            
            fprintf(fid,'Number of Material Property Mesh Resolutions: \n');
            fprintf(fid,'3 \n');
            fprintf(fid,'Number of material mesh structures \n');
            fprintf(fid,'3 \n');
            fprintf(fid,'Material mesh iteration stucture (one less than number of structures)  \n');
            fprintf(fid,['10 50 \n']);
            fprintf(fid,'Material Property Mesh Resolutions (x,y,z): \n');
            fprintf(fid,'%9.7f , %9.7f , %9.7f \n',vox);
            fprintf(fid,'%9.7f , %9.7f , %9.7f \n',2.*vox);
            fprintf(fid,'%9.7f , %9.7f , %9.7f \n',5.*vox);
            if(ifile==10)||(ifile==12) % Model 6, NITI aniso damping. Coarse mesh res at beginning to speed up convergence. 
                fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
                fprintf(fid,'3 3 3 3 3 \n');
                fprintf(fid,'3 3 3 3 3 \n');
                fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
                fprintf(fid,'2 2 2 3 3 \n');
                fprintf(fid,'2 2 2 3 3 \n');
                fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
                fprintf(fid,'1 1 1 3 3 \n');
                fprintf(fid,'1 1 1 3 3 \n');
                fprintf(fid,'Reconstruction Indicators: \n');
                fprintf(fid,'.true.,.true. \n');
                fprintf(fid,'.true.,.true. \n');
                fprintf(fid,'.true.,.true. \n');
                fprintf(fid,'.false.,.false. \n');
                fprintf(fid,'.false.,.false. \n');                
            elseif(ifile==11)||(ifile==13) % Model 6, NITI iso damping
                fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
                fprintf(fid,'3 3 3 3 3 \n');
                fprintf(fid,'3 3 3 3 3 \n');
                fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
                fprintf(fid,'2 2 2 3 3 \n');
                fprintf(fid,'2 3 3 3 3 \n');
                fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
                fprintf(fid,'1 1 1 3 3 \n');
                fprintf(fid,'1 3 3 3 3 \n');
                fprintf(fid,'Reconstruction Indicators: \n');
                fprintf(fid,'.true.,.true. \n');
                fprintf(fid,'.true.,.false. \n');
                fprintf(fid,'.true.,.false. \n');
                fprintf(fid,'.false.,.false. \n');
                fprintf(fid,'.false.,.false. \n');
            end
        end        
        fprintf(fid,'Property Estimate Variance Calculation: \n');
        fprintf(fid,'.false. \n');
        fprintf(fid,'Global forward solve to check displacement error \n');
        fprintf(fid,'.false. \n');
        fprintf(fid,'Zone Sizes (not including overlap factor, actual size = (1+2*ovlp)*siz) [x y z]: \n');
        fprintf(fid,'%12.5e,%12.5e,%12.5e',znedgelength);fprintf(fid,'\n');
        fprintf(fid,'Zone Overlap Percentage [x y z]: \n');
        fprintf(fid,'%4.3fd0,%4.3fd0,%4.3fd0',znovlp);fprintf(fid,'\n');
        fprintf(fid,'Iteration Limits [global zone]: \n');
        fprintf(fid,'%i,1 \n',globitr);
        fprintf(fid,'Minimum Parameter Update Size [global zone line]: \n');
        fprintf(fid,'%c','1.d-5, 1d-5, 1.d-5');fprintf(fid,'\n');
        fprintf(fid,'Number of Zone Iteration Structures: \n');
        fprintf(fid,'3 \n');
        if(ifile<10) % CG for isotropic models
            fprintf(fid,'Iteration Limits for Zone Iteration Structures (NOTE:  provide one less limit than # of structures!!!): \n');
            fprintf(fid,'10,150 \n');
            fprintf(fid,'Zone Iteration Structures [<# of CG iters> <# of GN iters> <# of QN iters> <# of line search iters>]: \n');
            fprintf(fid,'1,0,0,1 \n');
            fprintf(fid,'2,0,0,2 \n');
            fprintf(fid,'3,0,0,2 \n');
        else
            fprintf(fid,'Iteration Limits for Zone Iteration Structures (NOTE:  provide one less limit than # of structures!!!): \n');
            fprintf(fid,'10,70 \n');
            fprintf(fid,'Zone Iteration Structures [<# of CG iters> <# of GN iters> <# of QN iters> <# of line search iters>]: \n');
            fprintf(fid,'0,0,2,2 \n');
            fprintf(fid,'0,0,3,3 \n');
            fprintf(fid,'0,0,4,3 \n');
        end
        fprintf(fid,'Number of Processors per MUMPS Communicator \n');
        fprintf(fid,'1 \n');
        fprintf(fid,'Maximum Amount of RAM per Processor [MB] \n');
        fprintf(fid,'999999999 \n');
        fprintf(fid,'ooooooooooo  REGULARIZATION DESCRIPTORS  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n');
        fprintf(fid,'Regularization Indicators [<TV> <SF> <TK> <MQ> <JW> <constraints> <CG residual> <VH> <McG> <soft prior>]: \n');
        if(ifile==2)||(ifile==5)||(ifile==7)||(ifile==9)||(ifile==12)||(ifile==13)
            fprintf(fid,'.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true.,.true.,.true. \n');
        else
            fprintf(fid,'.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true.,.true.,.false. \n');
        end        
        fprintf(fid,'Number of constant regularization iterations (final N iterations use final regularization weights) \n');
        if(ifile<10)
            fprintf(fid,'%i \n',globitr-70);
        else
            fprintf(fid,'%i \n',40);
        end
        fprintf(fid,'Number of Parameters to Treat with Total Variation: \n');
        fprintf(fid,'%i \n',nprop);
        fprintf(fid,'TV Parameter List: \n');
        str=sprintf('%i,',1:nprop);fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'TV Delta Values: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('1.d-19,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('1.d-19,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'TV Initial Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('5.d-16,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('5.d-16,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'TV Final Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('5.d-16,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('5.d-16,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Number of Parameters to Treat with Spatial Filtering: \n');
        fprintf(fid,'%i \n',nprop);
        fprintf(fid,'SF Parameter List: \n');
        str=sprintf('%i,',1:nprop);fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'SF Gaussian filter initial Widths: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('0.003d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('0.003d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'SF Final Gaussian width: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('0.0015d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('0.0015d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Number of Parameters to Treat with Tikhonov Regularization: \n');
        fprintf(fid,'%i \n',nprop);
        fprintf(fid,'TK Parameter List: \n');
        str=sprintf('%i,',1:nprop);fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'TK Initial Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('1d-18,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('1d-18,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'TK Final Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('1d-18,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('1d-18,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Number of Parameters to Treat with Marquardt Regularization: \n');
        fprintf(fid,'%i \n',nprop);
        fprintf(fid,'Distance of alpha from 1 at which MQ weights are adjusted: \n');
        fprintf(fid,'0.25d0 \n');
        fprintf(fid,'MQ Parameter List: \n');
        str=sprintf('%i,',1:nprop);fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'MQ Weight Delta: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('0.5d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('0.5d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'MQ Initial Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('1d2,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('1d2,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'MQ Minimum Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('1d-11,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('1d-11,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Number of Parameters to Treat with Joachimowicz Regularization: \n');
        fprintf(fid,'%i \n',nprop);        
        fprintf(fid,'JW Parameter List: \n');
        str=sprintf('%i,',1:nprop);fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'JW Initial Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('5.d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('5.d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'JW Final Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('5.d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('5.d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Number of Parameters to Treat with Constraints: \n');
        fprintf(fid,'%i \n',nprop);        
        fprintf(fid,'Constraint Parameter List: \n');
        str=sprintf('%i,',1:nprop);fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'Constraint Weights <converges to exact solution as weight --> inf>: \n');
        str=sprintf('%c',repmat('1.d-14,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        if(ifile==1)||(ifile==2) % Iso incomp
            fprintf(fid,'Constraint Minimums: \n');
            fprintf(fid,'200.d0, 1.d0, 1000.d0, -5.d4, 1000.d0, 0.d0 \n');
            fprintf(fid,'Constraint Maximums: \n');
            fprintf(fid,'5.d5, 5.d5, 1000.d0, -1.d-8, 1.d12, 0.d0 \n');
        elseif(ifile==3)||(ifile==8)||(ifile==9) % iso comp
            fprintf(fid,'Constraint Minimums: \n');
            fprintf(fid,'100.d0, 1.d0, 1000.d0, -5.d4, 100.d0, 0.d0 \n');
            fprintf(fid,'Constraint Maximums: \n');
            fprintf(fid,'5.d5, 5.d5, 1000.d0, -1.d-8, 1.d12, 0.d0 \n');
        elseif(ifile==4)||(ifile==5)||(ifile==6)||(ifile==7) % Poro
            fprintf(fid,'Constraint Minimums: \n');
            fprintf(fid,'100.d0, 0.d0, 100.d0, 0.d0, -24.d0, -50.d0, 0.01d0, 0.d0 \n');
            fprintf(fid,'Constraint Maximums: \n');
            fprintf(fid,'5.d5, 5.d5, 5.d6, 5.d6, -5.d0, -5d0, 10.d0, 10.d0 \n');
        elseif(ifile==10)||(ifile==12)
            fprintf(fid,'Constraint Minimums: \n');
            fprintf(fid,'200.d0, 1.d0, 200.d0, 1.d0, 200.d0, 1.d0, 1000.d0, -5.d4, 1000.d0, 0.d0 \n');
            fprintf(fid,'Constraint Maximums: \n');
            fprintf(fid,'5.d5, 5.d5, 5.d5, 5.d5, 5.d5, 5.d5, 1000.d0, -1.d-8, 1.d12, 1.d12 \n');
        elseif(ifile==11)||(ifile==13)
            fprintf(fid,'Constraint Minimums: \n');
            fprintf(fid,'200.d0, 1.d0, -5.d0, -5.d0, -5.d0, -5.d0, 1000.d0, -5.d4, 1000.d0, 0.d0 \n');
            fprintf(fid,'Constraint Maximums: \n');
            fprintf(fid,'5.d5, 5.d5, 5.d0, 5.d0, 5.d0, 5.d0, 1000.d0, -1.d-8, 1.d12, 1.d12 \n');
        end
        if(ifile<10) % Isotropic models
            fprintf(fid,'CG/QN Residual Scaling Initial Weights: <1st line real, 2nd line imag> \n');
            str=sprintf('%c',repmat('0.04d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
            str=sprintf('%c',repmat('0.04d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
            fprintf(fid,'CG Residual Scaling Final Weights: <1st line real, 2nd line imag> \n');
            str=sprintf('%c',repmat('0.01d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
            str=sprintf('%c',repmat('0.01d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']); 
        elseif(ifile==10)||(ifile==12)
            fprintf(fid,'CG/QN Residual Scaling Initial Weights: <1st line real, 2nd line imag> \n');
            str=sprintf('%c',repmat('0.04d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
            str=sprintf('%c',repmat('0.04d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
            fprintf(fid,'CG Residual Scaling Final Weights: <1st line real, 2nd line imag> \n');
            str=sprintf('%c',repmat('0.01d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
            str=sprintf('%c',repmat('0.01d0,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']); 
        elseif(ifile==11)||(ifile==13) % Heavy weighting for phi and zeta updates. 
            fprintf(fid,'CG/QN Residual Scaling Initial Weights: <1st line real, 2nd line imag> \n');
            fprintf(fid,'0.04d0, 0.4d0, 0.4d0, 0.04d0, 0.04d0 \n');
            fprintf(fid,'0.04d0, 0.4d0, 0.4d0, 0.04d0, 0.04d0 \n');
            fprintf(fid,'CG Residual Scaling Final Weights: <1st line real, 2nd line imag> \n');
            fprintf(fid,'0.01d0, 0.1d0, 0.1d0, 0.01d0, 0.01d0 \n');
            fprintf(fid,'0.01d0, 0.1d0, 0.1d0, 0.01d0, 0.01d0 \n');
        end
        fprintf(fid,'Van Houten Regularization Level: \n');
        fprintf(fid,'1.2d0 \n');
        fprintf(fid,'Number of parameters to treat with soft prior: \n');
        fprintf(fid,'%i \n',nprop);        
        fprintf(fid,'SP Parameter List: \n');
        str=sprintf('%i,',1:nprop);fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'SP Initial Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('1.d-10,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('1.d-10,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'SP Final Weights: <1st line real, 2nd line imag> \n');
        str=sprintf('%c',repmat('1.d-10,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        str=sprintf('%c',repmat('1.d-10,',[1 nprop]));fprintf(fid,[str(1:end-1) '\n']);
        fprintf(fid,'SP start iteration delay (first N iterations do not use soft prior \n');
        fprintf(fid,'0 \n');
        fclose(fid);

    end

    