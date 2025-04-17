
function Plotpostprocess7p3(f,svfig)
% Prompt for file if none provided
if nargin == 0
    f = getfilename('Conv file to plot?', '*.prop.01.RE.convinfo');
end

% Trim path to base name
stm = f(1:end-20);

% Load reconstruction index
mshind = load([stm '.meshind']);
fid = fopen([stm '.reconind']);
rcnindtxt = fscanf(fid, '%c');
fclose(fid);
rcnindtxt = rcnindtxt(rcnindtxt == 'T' | rcnindtxt == 'F');
rcnind = rcnindtxt == 'T';

% Setup
np = 3;  % number of properties
nfig = 0;

% Loop through properties and real/imag components
for ip = 1:np
    ps = sprintf('%2.2i', ip);
    for jj = 1:2
        if jj == 1
            RI = 'RE'; RI_t = 'Real';
        else
            RI = 'IM'; RI_t = 'Imag';
        end

        if rcnind(2*(ip-1) + jj)
            fn = [stm '.prop.' ps '.' RI '.convinfo'];
            conv = load(fn);
            nitr = conv(end,1);

            nfig = nfig + 1;
            h(nfig) = figure;
            figinfo(nfig,1) = ip;
            figinfo(nfig,2) = jj;

            errorbar(1:nitr, conv(:,2), conv(:,3))
            hold on
            plot(1:nitr, conv(:,5), 'c.-', 1:nitr, conv(:,7), 'r.-')
            plot(1:nitr, conv(:,4), 'c.-', 1:nitr, conv(:,8), 'r.-')
            xlabel('Iteration'); ylabel('Value');
            title([RI_t ' - Property ' int2str(ip)]);
            legend('Mean','5%','95%','Min','Max');
        end
    end
end

%% Prompt to save figures
if nargin < 2
    svfig = input('Format to save figures <jpg, fig, tif, etc. | default = no save> >> ', 's');
end

if isempty(svfig)
    disp('No figures saved.')
else
    disp(['Saving figures as ' svfig ' images...'])
    for ii = 1:nfig
        if isgraphics(h(ii), 'figure')
            % Build file tag
            if figinfo(ii,2) == 1
                ftag = ['_conv_prop' sprintf('%2.2i', figinfo(ii,1)) '_Re'];
            else
                ftag = ['_conv_prop' sprintf('%2.2i', figinfo(ii,1)) '_Im'];
            end

            % Construct full file name
            savename = [stm ftag '.' svfig];

            % Save figure
            try
                saveas(h(ii), savename, svfig);
            catch ME
                warning('Failed to save figure %d: %s', ii, ME.message);
            end
        else
            warning('Skipping figure %d: invalid handle.', ii);
        end
    end
end

end



