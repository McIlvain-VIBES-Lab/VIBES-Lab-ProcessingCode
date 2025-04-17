% Plot a fluid flow field
% FluidFlowPlot('fluidsrcsim_tet4mm.hom.nod','fluidsrcsim_tet4mm.hom.elm','junk1.forward.flow','junk1.forward.dsp','junk1.forward.pre','../../NLI_MRE_Input.mat')
function FluidFlowPlot(nodf,elmf,flowf,dspf,presf,MREf)


junk=load(nodf);
x=junk(:,2);
y=junk(:,3);
z=junk(:,4);


junk=load(flowf);
qxr=junk(:,2);
qxi=junk(:,3);
qyr=junk(:,4);
qxi=junk(:,5);
qzr=junk(:,6);
qzi=junk(:,7);

load(MREf)
stack_dim=size(Ur(:,:,:,1));

flowstack=tet_output_general(flowf,nodf,elmf,2:7,stack_dim,voxsize_mm,'n');
dspstack=tet_output_general(dspf,nodf,elmf,2:7,stack_dim,voxsize_mm,'n');
prestack=tet_output_general(presf,nodf,elmf,2:3,stack_dim,voxsize_mm,'n');

montagestack(flowstack(1).vals);colorbar
title('flow x real')

montagestack(flowstack(3).vals);colorbar
title('flow y real')

montagestack(flowstack(5).vals);colorbar
title('flow z real')


montagestack(dspstack(1).vals);colorbar
title('disp x real')

montagestack(prestack(1).vals);colorbar
title('Pressure real')

k=1e-8;
eta=0.2;
rhoa=150;
rhof=1020;
omega=2*pi*1;

[dpxr,dpyr,dpzr]=gradient(prestack(1).vals);
[dpxi,dpyi,dpzi]=gradient(prestack(2).vals);

a=(k*1i*eta^2)/(1i*eta^2+omega*k*(rhoa+eta*rhof));
b=omega^2*rhof;

qxrg=a*(dpxr+b*dspstack(1).vals);
qxig=a*(dpxi+b*dspstack(2).vals);
qyrg=a*(dpyr+b*dspstack(3).vals);
qyig=a*(dpyi+b*dspstack(4).vals);
qzrg=a*(dpzr+b*dspstack(5).vals);
qzig=a*(dpzi+b*dspstack(6).vals);



montagestack(qxrg);colorbar
title('qxr: gradient based')
montagestack(qyrg);colorbar
title('qyr: gradient based')
montagestack(qzrg);colorbar
title('qzr: gradient based')

figure
quiver3(x,y,z,-qxr,-qyr,-qzr,2,'linewidth',2)

end

