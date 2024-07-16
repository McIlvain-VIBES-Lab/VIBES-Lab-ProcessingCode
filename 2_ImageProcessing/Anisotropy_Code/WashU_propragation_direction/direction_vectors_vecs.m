% script to plot weighted average propagation direction
% PVB Wash U 2013

% U3 - (NxMxP) field of complex coefficients (Re,Im) of fundamental
% component of scalar field (curl or displacement component)
% mask - mask array (NxMxP)

% For real data, get rid of NaNs from mask and data
% Also create mesh of cartesian coordinates
% Don't need to do this for artificial data
realdata = 1;
if realdata,
    load('Phantom_20130924_for_dirvec.mat');
    U3 = Wf;        % analyze z-component of displacment
    [N,M,P]=size(U3);
    [x,y,z]=meshgrid(1:M,1:N,1:P);
    nskip = 5;
else,
    % Pick the type of wave pattern
     imagetype = 'spheres';
    % imagetype = 'circles';
    % imagetype = 'ellipses';
%      imagetype = 'planewaves';
    
    % show movie if propmovie=1; don't if propmovie=0
    propmovie = 1;
    [U3,mask,x,y,z]=make_waves_propagating(imagetype, propmovie);
    nskip = 2;
end;

% get the direction vectors
[UUX,UUY,UUZ]=dirfilt_3D_vecs(U3,mask);

%% Make quiver plot in 3d
% usually don't plot all the vectors or it will take forever (skip some)

UX = UUX(1:nskip:end,1:nskip:end,1:nskip:end);
UY = UUY(1:nskip:end,1:nskip:end,1:nskip:end);
UZ = UUZ(1:nskip:end,1:nskip:end,1:nskip:end);

X = x(1:nskip:end,1:nskip:end,1:nskip:end);
Y = y(1:nskip:end,1:nskip:end,1:nskip:end);
Z = z(1:nskip:end,1:nskip:end,1:nskip:end);

figure(10)
quiver3(X,Y,Z,UX,UY,UZ)
axis equal
if realdata==0,
    title(['Ampl. weighted propag. dirs. for ',imagetype]);
else
   title('Ampl. weighted propag. dirs. for gel phantom');
end
% axis off