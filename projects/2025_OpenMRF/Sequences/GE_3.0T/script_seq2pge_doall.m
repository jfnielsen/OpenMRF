% Create .pge files and put them in OpenMRF_GE.tar
      
% -------------------------------------------------------------------------
% EDIT THIS SECTION AS NEEDED 
% -------------------------------------------------------------------------

isTest = false;  % if true, turn off PNS and max slew checks (for WTools)

% Local path containing all the .seq files you wish to convert to .pge
seqFilePath = '~/Downloads/OpenMRF/GE_OpenMRF_251117_maybe_final/';

% Output file name
tarFileName= 'OpenMRF-GE-' + replace(string(datetime), {':', ' '}, '-') + '.tar';

% Path on the scanner where the .pge files will reside, and the .entry
% file index corresponding to the first .pge file.
% These settings are used to generate the .entry files, which must be
% copied to /srv/nfs/psd/usr/psd/pulseq/v7/ on the scanner host computer.
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

% PNS channel/direction weights
PNSwt = (1-isTest)*[1 1 1];   

% -------------------------------------------------------------------------
% Convert .seq files to .pge format for the pge2 interpreter
% -------------------------------------------------------------------------

%pge2.utils.removeFiles('*.pge *.entry');

seqFilePath = pge2.utils.ensuretrailingslash(seqFilePath);
seqFilePath = pge2.utils.normalizepath(seqFilePath);

D = dir([seqFilePath '*.seq']);

% Initialize output .tar file
system('git rev-parse HEAD > commitID.txt');
system(sprintf('tar cf %s commitID.txt setup_4_seq2pge.m script_seq2pge_doall.m', tarFileName));
pge2.utils.removeFiles('commitID.txt');

for ii = 1:length(D)
    fn = replace(D(ii).name, '.seq', '');

    % Convert to Ceq sequence representation
    ceq = seq2ceq([seqFilePath fn '.seq'], 'usesRotationEvents', false);

    % Check PNS and b1/gradients against scanner limits,
    % and extract some sequence parameters.
    params = pge2.check(ceq, sysGE, 'wt', PNSwt);

    % Check accuracy of the ceq sequence representation against the .seq file
    sys = mr.opts('maxGrad', sysGE.g_max*10, 'gradUnit','mT/m', ...
              'maxSlew', sysGE.slew_max*10, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', sysGE.rf_dead_time, ...
              'rfRingdownTime', sysGE.rf_ringdown_time, ...
              'adcDeadTime', 10e-6, ...
              'adcRasterTime', sysGE.adc_raster_time, ...
              'rfRasterTime', sysGE.RF_UPDATE_TIME, ...
              'gradRasterTime', sysGE.GRAD_UPDATE_TIME, ...
              'blockDurationRaster', sysGE.GRAD_UPDATE_TIME);
    seq = mr.Sequence(sys);
    seq.read([seqFilePath fn '.seq']);
    xmlPath = []; % if nonempty, compare against WTools/Pulse Studio output
    pge2.validate(ceq, sysGE, seq, xmlPath, 'row', [], 'plot', false);

    % Write .pge file
    pislquant = 1;  % num ADC events for setting Rx gain in Auto Prescan
    pge2.writeceq(ceq, [fn '.pge'], 'pislquant', pislquant, 'params', ...
        pge2.utils.setfields(params, 'smax', (1-isTest) * params.smax));

    % Write .entry file
    entryFileNum = CV1 + ii - 1;
    pge2.writeentryfile(entryFileNum, fn, 'path', pgeFilePath);

    % Add files to tar file
    system(sprintf('tar --remove-files --append -f %s pge%d.entry %s', tarFileName, entryFileNum, [fn '.pge']));

    fprintf('\n\n\n%s\n', repmat('-', 1, 79));
end

