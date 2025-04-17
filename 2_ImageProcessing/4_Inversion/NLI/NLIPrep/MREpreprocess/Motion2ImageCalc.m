%Motion2Image
%Provided parms (MREMotionInfo structure) produces the correct motion
%transformation for : 
%-HEX & TET Mesh
%-Extrinsic and Intrinsic Motion Acquisition
% Scannertype==1 --> Philips
% Scannertype==2 --> SiemensProcessing and parms is DicomInfo structure
% built in SiemensProcessing.m



function [M2I]= Motion2ImageCalc(parms,scannertype)

if isempty(scannertype)
    scannertype=1; %Assumes Philips 
end

parms=parms(1); 

if scannertype==1;
    if parms.ExtInt==1 %Extrinsic Scan
        %% Computing Tpom matrix - going from magnetic to patient coordinate system
        % Head First
        if(strcmp(strcat(parms.patient_position),'Head First Prone')==1)
            tpo=[-1 0 0;0 -1 0;0 0 1];
            tpp=[0 1 0;-1 0 0;0 0 -1];
        end

        if(strcmp(strcat(parms.patient_position),'Head First Supine')==1)
            tpp=[0 1 0;-1 0 0;0 0 -1];
            tpo=[1 0 0;0 1 0;0 0 1];
        end

        if(strcmp(strcat(parms.patient_position),'Head First Rightcubitus')==1)
            tpp=[0 1 0;-1 0 0;0 0 -1];
            tpo=[0 -1 0;1 0 0;0 0 1] ;
        end

        if(strcmp(strcat(parms.patient_position),'Head First Leftcubitus')==1)
            tpp=[0 1 0;-1 0 0;0 0 -1];
            tpo=[0 1 0;-1 0 0;0 0 1]; 
        end

        % Feet First

        if(strcmp(strcat(parms.patient_position),'Feet First Prone')==1)
            tpo=[-1 0 0;0 -1 0;0 0 1];
            tpp=[0 -1 0;-1 0 0;0 0 1];
        end

        if(strcmp(strcat(parms.patient_position),'Feet First Supine')==1)
            tpp=[0 -1 0;-1 0 0;0 0 1];
            tpo=[1 0 0;0 1 0;0 0 1];
        end

        if(strcmp(strcat(parms.patient_position),'Feet First Rightcubitus')==1)
            tpp=[0 -1 0;-1 0 0;0 0 1];
            tpo=[0 -1 0;1 0 0;0 0 1];
        end

        if(strcmp(strcat(parms.patient_position),'Feet First Leftcubitus')==1)
            tpp=[0 -1 0;-1 0 0;0 0 1];
            tpo=[0 1 0;-1 0 0;0 0 1];
        end

        Tpom=tpo*tpp;

        %% computing the Tang matrix -- going from slice orientation to patient cord

        A1=parms.angulation;
        ap=A1(1);
        fh=A1(2);
        rl=A1(3);

        Trl=[1 0 0;0 cos(rl) -sin(rl);0 sin(rl) cos(rl)];
        Tap=[cos(ap) 0 sin(ap);0 1 0;-sin(ap) 0 cos(ap)];
        Tfh=[cos(fh) -sin(fh) 0;sin(fh) cos(fh) 0;0 0 1];

        Tang=Trl*Tap*Tfh;

        %% computing the  Tsom matrix -- going from image to slice orientation

        if(strcmp(parms.slice_orientation,'Transverse')==1)
            display('Transverse Scan')
            display(parms.preparation_dir)

        Tsom=[0 -1 0;-1 0 0;0 0 1];
        if(strcmp(parms.preparation_dir,'Right-Left')==1) 
            Tprep=[1 0 0; 0 1 0; 0 0 1];
            if isempty(parms.fatS)==1
            parms.fatS='P';
            end
        end


        if(strcmp(parms.preparation_dir,'Anterior-Posterior')==1)
             Tprep=[0 -1 0; 1 0 0; 0 0 1];

               if isempty(parms.fatS)==1
               parms.fatS='L';
               end 
        end
        display(parms.fatS)

        end

        if(strcmp(parms.slice_orientation,'Sagittal')==1)

            display('Sagital Scan')
            display(parms.preparation_dir)
        Tsom=[0 0 -1;0 -1 0;1 0 0];

        if(strcmp(parms.preparation_dir,'Feet-Head')==1)
             Tprep=[0 -1 0; 1 0 0; 0 0 1];

        end

        if(strcmp(parms.preparation_dir,'Anterior-Posterior')==1 )
             Tprep=[1 0 0; 0 1 0; 0 0 1];
        end


        display(parms.fatS);

        end

        %if(parms.slice_orientation==3)
        if (strcmp(parms.slice_orientation,'Coronal')==1);
            display('Coronal Scan')
            display(parms.preparation_dir)
            Tsom=[0 -1 0;0 0 1;1 0 0];
            if(strcmp(parms.preparation_dir,'Right-Left')==1)

               Tprep=[1 0 0; 0 1 0; 0 0 1];
               if isempty(parms.fatS)==1
               parms.fatS='F';
               end
            end

            if(strcmp(parms.preparation_dir,'Feet-Head')==1)
                Tprep=[0 -1 0; 1 0 0; 0 0 1];
                if isempty(parms.fatS)==1
                parms.fatS='L';
                end

            end

           display(parms.fatS) 

        end

        % Choice  m-filp, s-flip, s-flip based on slice orientation, fold over
        % direction and fat-shift direction

        % Transverse

        if(strcmp(parms.preparation_dir,'Right-Left')==1 && strcmp(parms.slice_orientation,'Transverse')==1)
            if(parms.fatS == 'R' || parms.fatS=='P')
            flips='p';
            end

            if(parms.fatS =='L' || parms.fatS=='F' || parms.fatS=='A')
            flips='m';
            end

            if(parms.fatS =='H')
            flips='s';
            end
        end

        if(strcmp(parms.preparation_dir,'Anterior-Posterior')==1 && strcmp(parms.slice_orientation,'Transverse')==1)    
            if(parms.fatS=='R' || parms.fatS=='F' || parms.fatS=='A')
              flips='m';
            end

            if(parms.fatS=='L' || parms.fatS=='P')
              flips='p';
            end

            if(parms.fatS=='H')
              flips='s';      
            end
        end

        % Sagital

        if(strcmp(parms.preparation_dir,'Feet-Head')==1 && strcmp(parms.slice_orientation,'Sagittal')==1)

          if(parms.fatS=='L' || parms.fatS=='H' || parms.fatS=='A')
           flips='m';
          end

          if(parms.fatS=='F' || parms.fatS=='P')
           flips='p';
          end

          if(parms.fatS=='R')
          flips='s';
          end

        end  
        if(strcmp(parms.preparation_dir,'Anterior-Posterior')==1 && strcmp(parms.slice_orientation,'Sagittal')==1)    

          if(parms.fatS=='L' || parms.fatS=='H' || parms.fatS=='P')
           flips='m';
          end

          if(parms.fatS=='A' || parms.fatS=='F')
           flips='p';
          end

          if(parms.fatS=='R')
          flips='s';
          end

        end

        % Coronal
        if(strcmp(parms.preparation_dir,'Feet-Head')==1 && strcmp(parms.slice_orientation,'Coronal')==1)

          if(parms.fatS=='R' || parms.fatS=='H' || parms.fatS=='A')
           flips='m';
          end

          if(parms.fatS=='L' || parms.fatS=='F')
           flips='p';
          end

          if(parms.fatS=='P')
          flips='s';
          end 

        end
        if(strcmp(parms.preparation_dir,'Right-Left')==1 && strcmp(parms.slice_orientation,'Coronal')==1)

          if(parms.fatS=='L' || parms.fatS=='H' || parms.fatS=='A')
           flips='m';
          end

          if(parms.fatS=='R' || parms.fatS=='F')
           flips='p';
          end

          if(parms.fatS=='P')
          flips='s';    
          end

        end


        if (flips=='m')   
        Tfsd=[-1 0 0;0 1 0;0 0 1];
        end

        if(flips=='p')
        Tfsd=[1 0 0;0 -1 0;0 0 1];
        end

        if(flips=='s')
        Tfsd=[1 0 0;0 1 0;0 0 -1];
        end

        Txyz=Tprep*Tfsd;

        LHS=inv(Tpom)*Tang*Tsom;
        RHS=inv(Tpom)*Tang*Tsom*Txyz;
        NWV2ijk=[0 -1 0; -1 0 0; 0 0 1];

        M2I = NWV2ijk*inv(LHS)*RHS;

    elseif parms.ExtInt==2 %Intrinsic 2DQFLOW / 4DQFLOW

        if strcmp(parms(1).slice_orientation,'Transverse')
            M2I=[0	1	0; 1	0	0; 0	0	1]; 
        elseif strcmp(parms(1).slice_orientation,'Coronal')
            M2I=[0	1	0; 0	0	1;-1	0	0];
        elseif strcmp(parms(1).slice_orientation,'Sagittal')
            M2I=[0	0	-1;  0	1	0;-1	0	0];
        end
    end
elseif scannertype==2 
    if strcmp(parms.Private_0051_100e,'Tra') && strcmp(parms.InPlanePhaseEncodingDirection,'COL'); 
        M2I=[0 1 0; 1 0 0; 0 0 1];
    elseif strcmp(parms.Private_0051_100e,'Tra') && strcmp(parms.InPlanePhaseEncodingDirection,'ROW'); 
        M2I=[0 -1 0; -1 0 0; 0 0 1];
    elseif strcmp(parms.Private_0051_100e,'Cor') && strcmp(parms.InPlanePhaseEncodingDirection,'COL'); 
        M2I=[0 0 1; 1 0 0; 0 1 0];
    elseif strcmp(parms.Private_0051_100e,'Cor') && strcmp(parms.InPlanePhaseEncodingDirection,'ROW'); 
        M2I=[0 0 -1; -1 0 0; 0 1 0];        
    elseif strcmp(parms.Private_0051_100e,'Sag') && strcmp(parms.InPlanePhaseEncodingDirection,'COL'); 
        M2I=[0 0 1; 0 -1 0; 1 0 0];
    elseif strcmp(parms.Private_0051_100e,'Sag') && strcmp(parms.InPlanePhaseEncodingDirection,'ROW'); 
        M2I=[0 0 -1; 0 1 0; 1 0 0];    
    end
end

%Tag of 1 or 0 for 

