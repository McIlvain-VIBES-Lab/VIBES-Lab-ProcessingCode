function [U3,mask,x,y,z]=make_waves_propagating(imagetype, propmovie)
% make waves in 3D
% test LFE 3d
% show propagating waves
% PVB Wash U 2014

% Pick the type of wave pattern
% imagetype = 'spheres';
% imagetype = 'circles';
% imagetype = 'ellipses';
% imagetype = 'planewaves';

% show movie if 1; don't of 0
% propmovie = 1;

% size of artificial domain
Lx = 1;
Ly = 1;
Lz = 1;

% discretization parameters
Nx = 32;
Ny = Nx;
Nz = Nx;

dx = Lx/Nx;
dy = dx;
dz = dx;
[x,y,z]=meshgrid(-Lx/2:dx:(Lx/2-dx),-Ly/2:dy:(Ly/2-dy),-Lz/2:dz:(Lz/2-dz));

switch imagetype,
    
    case 'spheres',
        % images
        R=sqrt(x.^2+y.^2+z.^2);
        kr = 5;  % wavenumber
        U3=exp(1i*(kr*2*pi*R));
        mask = R<Lx/2;
        U3=U3.*mask;
    case 'circles',
        % images
        R=sqrt(x.^2+y.^2);
        kr = 5;  % wavenumber
        U3=(exp(1i*(kr*2*pi*R)));
        mask = R<Lx/2;
        U3=U3.*mask;
    case 'ellipses',
        % images
        ka = 3;  % wavenumber
        kb = 5;  % wavenumber
        R=sqrt((ka*x).^2+(kb*y).^2);
        % U3=(exp(1i*(2*pi*R)));        %inward
        U3=conj(exp(1i*(2*pi*R)));      % outward
        mask = R<2*Lx;
        U3=U3.*mask;
    case 'planewaves',
        kx = -4;
        ky = -2;
        kz = 2;
        U3=exp(1i*(kx*2*pi*x + ky*2*pi*y + kz*2*pi*z));
        mask = abs(x)<Lx/2;
        U3=U3.*mask;
end;

figure(1)
for n=1:32,
    subplot(6,6,n),
    imagesc(real(U3(:,:,n)),[-1 1])
    axis equal, axis off, axis xy
end;


if propmovie,       % show two cycles of wave propagation
    figure(2)
    t=pi/16:pi/16:2*pi;
    for n=1:32,
        u(:,:,:,n)=U3*exp(2*1i*t(n));
    end;
    for n=1:32,
        imagesc(real(u(:,:,16,n)),[-1 1]);
        axis equal,
        axis off
        axis xy
        mov(:,n)=getframe;
    end;
end;

% 3D FFT of volume, then take a slice of FFT
f2n = fftn(U3)/(32^3);
f2ns = fftshift(f2n);
f2s = squeeze(f2ns(:,:,17));
axf = [-16:15]/16;

figure(3),
switch imagetype,
    case 'planewaves',
        subplot(2,2,1),imagesc(real(U3(:,:,16)),[-1 1]), axis equal, axis off, axis xy; title(imagetype);
        subplot(2,2,3),imagesc(abs(f2s),[0 1]), colormap gray, axis equal, axis off, axis xy; title(['FFT of ',imagetype]);
    case 'spheres',
        subplot(2,2,2),imagesc(real(U3(:,:,16)),[-1 1]), axis equal, axis off,  axis xy; title(imagetype);
        subplot(2,2,4),imagesc(abs(f2s),[0 0.25]), colormap gray, axis equal, axis off, axis xy; title(['FFT of ',imagetype]);
    case 'circles',
        subplot(2,2,2),imagesc(real(U3(:,:,16)),[-1 1]), axis equal, axis off, axis xy; title(imagetype);
        subplot(2,2,4),imagesc(abs(f2s),[0 0.25]), colormap gray, axis equal, axis off, axis xy;  title(['FFT of ',imagetype]);
    case 'ellipses',
        subplot(2,2,2),imagesc(real(U3(:,:,16)),[-1 1]), axis equal, axis off, axis xy; title(imagetype)
        subplot(2,2,4),imagesc(abs(f2s),[0 0.25]), colormap gray, axis equal, axis off, axis xy;  title(['FFT of ',imagetype]);
end;






