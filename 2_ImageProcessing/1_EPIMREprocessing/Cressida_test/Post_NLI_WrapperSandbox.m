%Post-NLI processing wrapper script
%Sept 17th 2025, Cressida Michaloski and Teah Serani   

%Section 1: Check for NLI outputs
    CheckForNLIOutputs();

%Section 2: Register files on Insomnia
    system(['bash /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/Cressida_test/reg_insomniasandbox.sh'])
    %add insomnia script for registration somehow if not already in bash script?

%Section 3: Pull freesurfer
