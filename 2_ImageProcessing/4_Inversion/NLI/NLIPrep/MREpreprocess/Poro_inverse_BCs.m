% Scriptfile Poro_inverse_BCs.m
% ************************************************************************
% This file reads in the Boundary Nodes-file of a mesh and adds boundary 
% conditions to the nodes selected by an
% appropriate condition.
%
% Input: Converted files from mesh:
%        - node file (*.nod)
%        - bndfile (*.bnd)
%
% Output: A file 'name.bcs' containing the boundary condition farmatted:
% <BC-#> <Node-#> <x-type> <jnk> <xvalue>  <y-type> <jnk> <yvalue> <z-type> <jnk> <zvalue>
%
% Original Author: Ashton Peters
% Modified: Hans-Uwe Berger
%   14.07.2005: 
%     - Introduced the concideration of boundary nodes only
%     - Integrated 3D view of Bnodes enabled to be turned and B-values
%       possibly changed when accidently chosen wrong.
%   26.08.2005: 
%     - Generalized in such, that program works on either tetrahedral or 
%       hexahedral elements, as long as the format of the mesh source files 
%       is given as: 
%       - bel-file: <Bel-#> <Bnode-# 1> <Bnode-# 2> ... <Bnode-# n> <elm-#>
%       - nod-file:  <Node-#> <x-pos.> <y-pos.> <z-pos.>
% Modified 5/18/2017
%     - Only does pressure BCs for inverse problem
%     
% ************************************************************************
clear all; close all; clc;



% 1.) Get mesh-data:
% ==================
% Read in mesh-file stem:
disp(' ');
disp(' Script-File addbelBConds.m');
disp(' ');
disp(' This script reads in the Boundary Nodes-file of a mesh');
disp(' and adds boundary conditions to selected nodes.');
disp(' ');
disp(' Files:');
eval('ls *.nod')
disp(' ');
inp = input([' Please enter the mesh-file stem: >> '],'s');
if isempty(inp)
    disp(' ERROR - no such file, please do again: ')
    inp = input([' File-stem: >> '],'s');
end
disp(' ');

nod = load([inp,'.nod']);
bel = load([inp,'.bnod']);

disp(' + .nod-file loaded');
disp(' + .bnd-file loaded');
disp(' ');




Bnodes = (bel(:,2));
clear nodes Node bel



% 3.) Read out position data and plot:
% ====================================
%nn=size(nod,1);
x = nod(Bnodes,2);
y = nod(Bnodes,3);
z = nod(Bnodes,4);
clear nod

figure(1)
set(1,'WindowStyle','docked');
plot3(x,y,z,'.','MarkerEdgeColor',[0.8 0.8 0.8]);
axis equal; view([1.5 1 1]); 
hold on;
xlabel('x'); ylabel('y'); zlabel('z');



% 4.) Loop to create BCs:
% =======================
cont=1;
nbcs=0;
nbcsloc=0;
while cont==1
    bfind = 0;
    disp(' To find appropriate boundaries, enter a condition of the form "x > 0.039",');
    disp(' "z < 0.01", or similar.');
    while bfind == 0
        str = input(' Condition: >> ','s');
        eval(['bnodes = find(',str,');'])
        disp(' ');
        
        % Update the plot!
        plot3(x(bnodes),y(bnodes),z(bnodes),'y.');
        
        %figure(1)
        if exist('specify') == 0
            specify = 'y';
        end
        inp = input([' Take these Boundaries to specify conditions ([y]/[n], default: ',specify,' )? >> '],'s');
        disp(' ');
        if isempty(inp);
            inp = specify;
        end
        specify = inp;
        if inp == 'y'
            bfind = 1;
        else
            plot3(x(bnodes),y(bnodes),z(bnodes),'.','MarkerEdgeColor',[0.8 0.8 0.8]);
        end
    end
    
    % Specify typ:
    if exist('typ') == 0
        typ = [1];
    end    
    disp(' Please enter the Type of BC for the nodes specified in the form');
    disp(' "[xtyp ytyp ztyp]", where Type "1" = fixed pressure, Type "2" = Flow,');
    inp = input(' Default is [1]: >> ');
    disp(' ');
    if isempty(inp)
        inp = typ;
    end
    typ = inp;
    
    % Specify Real part of Value
    if exist('Revalue') == 0
        Revalue = [0];
    end
    disp(' Please enter the real part of the BC-Value for the nodes specified in the form');
    inp = input([' "[xval yval zval pval]", Default is [',num2str(Revalue),'] : >> ']);
    disp(' ');
    if isempty(inp)
        inp = Revalue;
    end
    Revalue = inp;
    
    % Specify Real part of Value
    if exist('Imvalue') == 0
        Imvalue = [0];
    end
    disp(' Please enter the Imaginary part of the BC-Value for the nodes specified in the form');
    inp = input([' "[xval yval zval pval]", Default is [',num2str(Imvalue),'] : >> ']);
    disp(' ');
    if isempty(inp)
        inp = Imvalue;
    end
    Imvalue = inp;
    
    % Update Plot:
    if typ==1
        plot3(x(bnodes),y(bnodes),z(bnodes),'r.');
    elseif typ==2 
        plot3(x(bnodes),y(bnodes),z(bnodes),'g.');
    else
        plot3(x(bnodes),y(bnodes),z(bnodes),'b.');
    end
    disp(' Plot updated.');
    disp(' ');
    
    % Save BC's as a variable bcs:
    for i=(nbcs+1):(nbcs+size(bnodes))
        bcs(i,1)=i;
        bcs(i,2)=Bnodes(bnodes(i-nbcsloc));
        bcs(i,3:5)=[typ Revalue Imvalue];
        nbcs=nbcs+1;
    end
    nbcsloc=nbcs;

    % Add another BC -> Update Loop!
    if exist('another') == 0
        another = 'y';
    end
    inp = input([' Want to add another BC ? ([y]/[n], Default: ',another,') >> '],'s');
    if isempty(inp)
        inp = another;
    end
    another = inp;

    if another == 'y'
        cont = 1;
    else
        cont = 0;
    end
end

% 4.2) Delete Duplicates !!!
% ==========================


% 5.) Write specified Bondary Conditions to file:
% ===============================================
if exist('outfile') == 0
    outfile = 'output.bcs';
end
disp(' Please enter a name for the output file to be written,');
inp = input([' Default is ',outfile,': >>  '],'s');
if isempty(inp)
    inp = outfile;
end
outfile = inp;

fid = fopen(outfile,'w');
fprintf(fid,'%7i  %7i  %i  %f  %f \n',bcs');
fclose(fid);

disp(' ');
disp([' File "',outfile,'" successfully written!']);
fprintf(' Total boundary conditions: %7.0f\n',nbcs);
disp(' ');
disp(' ***** END OF PROGRAM *****');




