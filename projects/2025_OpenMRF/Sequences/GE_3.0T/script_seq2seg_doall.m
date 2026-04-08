% Create .pge files and put them in OpenMRF_GE.tar
      
% -------------------------------------------------------------------------
% EDIT THIS SECTION AS NEEDED 
% -------------------------------------------------------------------------

isTest = false;  % if true, turn off PNS and max slew checks (for WTools)

% Local path containing all the .seq files you wish to convert to .pge
seqFilePath = '~/Downloads/OpenMRF/';

% Output file name
tarFileName= 'OpenMRF-GE-' + replace(string(datetime), {':', ' '}, '-') + '.tar';

% Path on the scanner where the .pge files will reside, and the .entry
% file index corresponding to the first .pge file.
% These settings are used to generate the .entry files, which must be
% copied to /srv/nfs/psd/usr/psd/pulseq/v7/ on the scanner host computer.
CV1 = 721;               % any non-negative integer
pgeFilePath = '/srv/nfs/psd/usr/psd/pulseq/v7/sequences/OpenMRF';

% Scanner hardware settings
psd_rf_wait = 100e-6;    % RF-gradient delay, scanner specific (s)
psd_grd_wait = 100e-6;   % ADC-gradient delay, scanner specific (s)
b1_max = 0.25;           % Gauss
g_max = 5;               % Gauss/cm
slew_max = 20;           % Gauss/cm/ms
coil = 'xrm';            % MR750. See pge2.opts()
sys_ge = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

% PNS channel/direction weights
PNSwt = (1-isTest)*[1 1 1];   

% -------------------------------------------------------------------------
% Convert .seq files to .pge format for the pge2 interpreter
% -------------------------------------------------------------------------

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
    psq = pulseg.fromSeq([D(ii).folder '/' seq_name '.seq']);   % ,'usesRotationEvents', false);

    % Check PNS and b1/gradients against scanner limits,
    % and extract some sequence parameters.
    PNSwt = 0.95*[1 1 1];   % directional PNS weights, see pge2.pns()
    params = pge2.check(ceq, sys_ge, 'PNSwt', PNSwt);

    % Check accuracy of the ceq sequence representation against the .seq file
    sys = mr.opts('maxGrad', sys_ge.g_max*10, 'gradUnit','mT/m', ...
              'maxSlew', sys_ge.slew_max*10, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', sys_ge.rf_dead_time, ...
              'rfRingdownTime', sys_ge.rf_ringdown_time, ...
              'adcDeadTime', sys_ge.adc_dead_time, ...
              'adcRasterTime', sys_ge.adc_raster_time, ...
              'rfRasterTime', sys_ge.RF_UPDATE_TIME, ...
              'gradRasterTime', sys_ge.GRAD_UPDATE_TIME, ...
              'blockDurationRaster', sys_ge.GRAD_UPDATE_TIME);
    seq = mr.Sequence(sys);
    seq.read([seqFilePath fn '.seq']);
    xmlPath = []; % if nonempty, compare against WTools/Pulse Studio output
    pge2.validate(ceq, sys_ge, seq, xmlPath, 'row', [], 'plot', false);

    % save to .mat file
    save([seq_name '.mat'], 'psq', 'params');
    system(sprintf('tar --append -f %s %s.mat', tarfn, seq_name));

    % Write .pge file
    pislquant = 1;  % num ADC events for setting Rx gain in Auto Prescan
    pge2.writeceq(ceq, [fn '.pge'], 'pislquant', pislquant, 'params', ...
        pge2.utils.setfields(params, 'smax', (1-isTest) * params.smax));

    % Write .entry file
    entryFileNum = CV1 + ii - 1;
    pge2.writeentryfile(entryFileNum, fn, 'path', pgeFilePath);
    system(sprintf('tar --append -f %s pge%d.entry', tarFileName, CV1+ii-1));

    % Add files to tar file
    system(sprintf('tar --remove-files --append -f %s pge%d.entry %s', tarFileName, entryFileNum, [fn '.pge']));

    fprintf('\n\n\n%s\n', repmat('-', 1, 79));
end

