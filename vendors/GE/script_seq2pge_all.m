% Create .pge files and put them in OpenMRF.tar
      
% -------------------------------------------------------------------------
% EDIT THIS SECTION AS NEEDED 
% -------------------------------------------------------------------------
% Local path containing all the .seq files you wish to convert to .pge
seqFilePath = '~/Downloads/OpenMRF/GE_OpenMRF_251117_maybe_final/';

% Output file name
tarFileName = 'OpenMRF_GE.tar'; 

% Path on the scanner where .pge files will reside,
% and .entry file number corresponding to the first .pge file.
% These settings are used to create the .entry files, that you must copy
% to /srv/nfs/psd/usr/psd/pulseq/v7/ on the scanner host computer.
CV1 = 721;    
pgeFilePath = '/srv/nfs/psd/usr/psd/pulseq/v7/sequences/OpenMRF';

% Scanner hardware settings
psd_rf_wait = 100e-6;    % RF-gradient delay, scanner specific (s)
psd_grd_wait = 100e-6;   % ADC-gradient delay, scanner specific (s)
b1_max = 0.25;           % Gauss
g_max = 5;               % Gauss/cm
slew_max = 20;           % Gauss/cm/ms
coil = 'xrm';            % MR750. See pge2.opts()
sysGE = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

% -------------------------------------------------------------------------
% Convert .seq files to .pge format for the pge2 interpreter
% -------------------------------------------------------------------------

seqFilePath = pge2.utils.ensuretrailingslash(seqFilePath);
seqFilePath = pge2.utils.normalizepath(seqFilePath);

D = dir([seqFilePath '*.seq']);

% Initialize output .tar file
system(sprintf('rm -f %s', tarFileName));
system(sprintf('tar cf %s main.m', tarFileName));

for ii = 1:length(D)
    fn = replace(D(ii).name, '.seq', '');

    if contains(fn, 'wasabi'), continue, end

    % Convert to Ceq sequence representation
    ceq = seq2ceq([seqFilePath fn '.seq'], 'usesRotationEvents', false);

    % Check PNS and b1/gradients against scanner limits,
    % and extract some sequence parameters.
    PNSwt = 0*[1 1 1];   % PNS channel/direction weights
    params = pge2.check(ceq, sysGE, 'wt', PNSwt);
    %params.smax = 1;   % for simulating in WTools

    % Check accuracy of the ceq sequence representation against the .seq file
    seq = mr.Sequence();
    seq.read([seqFilePath fn '.seq']);
    warning('OFF', 'mr:restoreShape');  % turn off Pulseq warning for spirals
    xmlPath = []; % if nonempty, compare against WTools/Pulse Studio output
    pge2.validate(ceq, sysGE, seq, xmlPath, 'row', [], 'plot', false);

    % Write .pge file
    pislquant = 1;  % num ADC events for setting Rx gain in Auto Prescan
    pge2.writeceq(ceq, [fn '.pge'], 'pislquant', pislquant, 'params', params);

    % Write .entry file
    entryFileNum = CV1 + ii - 1;
    pge2.writeentryfile(entryFileNum, fn);

    % Add files to tar file
    system(sprintf('tar --append -f %s pge%d.entry %s', tarFileName, entryFileNum, [fn '.pge']));
end
