function [Gpnan_all,Gdpnan_all,nrenan_all] = LDI_ShearModulus_Plot(LDIsummaryfile,cm)
%LDI_ShearModulus_Plot : makes a plot of all the data in the summary file
% inputs: 
%    summaryfile: output file from run_LDI_fun.m
%    cm: colormap (if not specified, cm = UIUCmap)
% written by RJO 11/21/16

load(LDIsummaryfile)
[~,nm,nmext]=fileparts(LDIsummaryfile);
if ~strcmp(nmext,'.mat')
    error('wrong file type')
end
if nargin < 2,
    cm=UIUCelastogram;
else
    if size(cm,2)~=3;
        error('invalid colormap')
    end
end
if length(fitrange) > 1, 
    use=1:length(fitrange);
    fitn=2*fitrange+1;
    dim4label=['fit kern: ',sprintf('%d /',fitn)];
    dim4label=dim4label(1:end-1); %drop final slash
elseif length(Delta)> 1,
    use=1:length(Delta)
    dim4label=['k ,',sprintf('%d/',Delta)];
    dim4label=dim4label(1:end-1); %drop final slash
else
    use=1;
    dim4label='n/a';
end
%if ~exist('Gpnan_all','var'), 
for k=use,
    load(savefilename{k})
    Gpnan=Gp;Gpnan(Gp==0)=NaN;
    Gdpnan=Gdp;Gdpnan(Gdp==0)=NaN;
    nrenan=nre;nrenan(nre==0)=NaN;
    Gpnan_all(:,:,:,k)=Gpnan; Gdpnan_all(:,:,:,k)=Gdpnan; nrenan_all(:,:,:,k)=nrenan;
end
save(LDIsummaryfile,'-append','Gpnan_all','Gdpnan_all','nrenan_all')
%end % 
disprange=[0 max(abs(1e-3*Gpnan_all(:)))];
vis5d_profile(cat(5,1e-3*Gpnan_all,1e-3*Gdpnan_all,nrenan_all),colormap(cm),nm,disprange,{dim4label,'Re/Im/nre'})
end

