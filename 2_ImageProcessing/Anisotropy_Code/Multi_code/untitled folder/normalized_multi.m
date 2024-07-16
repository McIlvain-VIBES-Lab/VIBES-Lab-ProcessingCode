function [MuN] = normalized_multi(xy,z)
set(0,'DefaultFigureWindowStyle','docked') %dock all figures
load(sprintf('/Volumes/CLJ-001/Group/drsmitty/Code/Matlab/mre_colormaps.mat'))

addpath(sprintf('%s/AP/PostInv',pwd)); 
load('MuMasked.mat'); MuAP = MuNew; 
rmpath(sprintf('%s/AP/PostInv',pwd))

addpath(sprintf('%s/LR/PostInv',pwd)); 
load('MuMasked.mat'); MuLR = MuNew; 
rmpath(sprintf('%s/LR/PostInv',pwd))

addpath(sprintf('%s/SI/PostInv',pwd)); 
load('MuMasked.mat'); MuSI = MuNew; 
rmpath(sprintf('%s/SI/PostInv',pwd))

MuN = zeros(size(MuAP));

for i = 1:xy
    for j = 1:xy
        for k = 1:z
            MuN(i,j,k) = min([MuAP(i,j,k) MuLR(i,j,k) MuSI(i,j,k)]);
        end
    end
end

cd('./AP/PostInv')
save BaseMu.mat MuN
MuNorm = MuAP./MuN;
%montagestack(MuNorm); caxis([0 2.5]); set(gcf,'Colormap',stiff_color)
MuNorm = silencer(MuNorm,xy,z);
save MuNorm.mat MuNorm
montagestack(MuNorm); caxis([0 2]); set(gcf,'Colormap',stiff_color)
cd ..; cd ..

cd('./LR/PostInv')
save BaseMu.mat MuN
MuNorm = MuLR./MuN;
%montagestack(MuNorm); caxis([0 2.5]); set(gcf,'Colormap',stiff_color)
MuNorm = silencer(MuNorm,xy,z);
save MuNorm.mat MuNorm
montagestack(MuNorm); caxis([0 2]); set(gcf,'Colormap',stiff_color)
cd ..; cd ..

cd('./SI/PostInv')
save BaseMu.mat MuN
MuNorm = MuSI./MuN;
%montagestack(MuNorm); caxis([0 2.5]); set(gcf,'Colormap',stiff_color)
MuNorm = silencer(MuNorm,xy,z);
save MuNorm.mat MuNorm
montagestack(MuNorm); caxis([0 2]); set(gcf,'Colormap',stiff_color)
cd ..; cd ..

montagestack(MuN); caxis([0 5000]); set(gcf,'Colormap',stiff_color)

%%
function [silenced] = silencer(Matrix,xy,z)
silenced = zeros(size(Matrix));
for i = 1:xy
    for j = 1:xy
        for k = 1:z
            %Matrix(i,j,k)
            if Matrix(i,j,k) > 1.2
                silenced(i,j,k) = Matrix(i,j,k);
            else 
                silenced(i,j,k) = 0;
            end
        end
    end
end
%% 