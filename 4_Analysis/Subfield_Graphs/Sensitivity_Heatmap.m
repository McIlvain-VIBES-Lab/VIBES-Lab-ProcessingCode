%% Heatmap for Sensitivity Effect Sizes

figure;
Value = [0.103 0.269; 0.078 0.216] 

finalValue = Value; 
% Frequency = columns; Resolution = rows ( 1.5, 2, 2.5, 3mm)
xvalues = {'1.5', '0.9'};
yvalues = {'12', '10'};
h = heatmap(xvalues,yvalues,finalValue);

h.Title = 'Effect Size Values - Stiffness';
h.XLabel = 'Spatial Filtering Value';
h.YLabel = 'Soft Prior Regulartization Weighting';
h.MissingDataLabel = 'No data';
h.MissingDataColor = [0 1 1];
%h.Colormap = winter;
h.GridVisible = 'on'; % removes axis lines
h.CellLabelColor = 'none'; % removes data labels
set(gca, 'fontsize',14)
caxis(h, [0.05 0.3]);
%print('-dpng','180410_Contrast_Heatmap_slice24')


figure;
finalValue2 = [0.625 0.543; 0.63 0.546];
xvalues2 = {'1.5', '0.9'};
yvalues2 = {'12', '10'};
h = heatmap(xvalues2,yvalues2,finalValue2);
h.Title = 'Effect Size Values - Damping Ratio';
h.XLabel = 'Spatial Filtering Value';
h.YLabel = 'Soft Prior Regulartization Weighting';
h.MissingDataLabel = 'No data';
h.MissingDataColor = [0.9 0.9 0.9];
%h.Colormap = winter;
h.GridVisible = 'on'; % removes axis lines
h.CellLabelColor = 'none'; % removes data labels
set(gca, 'fontsize',14)
caxis(h, [0.5 0.65]);


