% Create .pge files and put them in OpenMRF_GE.tar
      
% -------------------------------------------------------------------------
% EDIT THIS SECTION AS NEEDED 
% -------------------------------------------------------------------------

is_test = true;  % if true, turn off PNS and max slew checks (for WTools)
is_traj = true;

% Local path containing all the .seq files you wish to convert to .pge
%seq_file_path = '~/Downloads/OpenMRF/mrf/';
if is_traj
    seq_file_path = '~/Downloads/OpenMRF/traj_ge/';
    tar_file_name= 'OpenMRF-traj_GE-' + replace(string(datetime), {':', ' '}, '-') + '.tar';
    opuser1 = 731;   % .entry file number for first scan (any positive integer)
else
    seq_file_path = '~/Downloads/OpenMRF/mrf/';
    tar_file_name= 'OpenMRF-GE-' + replace(string(datetime), {':', ' '}, '-') + '.tar';
    opuser1 = 721;
end

% Scanner hardware settings (for checking compatibility)
psd_rf_wait = 100e-6;    % RF-gradient delay, scanner specific (s)
psd_grd_wait = 100e-6;   % ADC-gradient delay, scanner specific (s)
b1_max = 0.25;           % Gauss
g_max = 5;               % Gauss/cm
slew_max = 20;           % Gauss/cm/ms
coil = 'xrm';            % MR750. See pge2.opts()
sys_ge = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

% PNS channel/direction weights
PNSwt = 0.9 * (1-is_test) * [1 1 1];   

% Pulseq scan list file (for interactive FOV prescription on scanner)
scans_list_file = 'pulseq_scans.list';

% -------------------------------------------------------------------------
% Convert each .seq file to PulSeg representation and save as .mat file
% -------------------------------------------------------------------------

seq_file_path = ensuretrailingslash(seq_file_path);
seq_file_path = normalizepath(seq_file_path);

D = dir([seq_file_path '*.seq']);

% Initialize output .tar file
system('git rev-parse HEAD > commitID.txt');
system(sprintf('tar cf %s commitID.txt setup_4_seq2pge.m script_seq2seg_doall.m', tar_file_name));
removefiles('commitID.txt');

% Initialize pulseq_scans.list file
fid = fopen(scans_list_file, 'w');
fprintf(fid, '# opuser1\tscan\n');   

for ii = 1:length(D)
    seq_name = replace(D(ii).name, '.seq', '');

    % Update .list file
    opuser1 = opuser1 + (ii > 1);
    fprintf(fid, '%d\t%s.mat\n', opuser1, seq_name);

    % Convert to PulSeg sequence representation
    psq = pulseg.fromSeq([D(ii).folder '/' seq_name '.seq']);   % ,'usesRotationEvents', false);

    % Check PNS and b1/gradients against scanner limits,
    % and extract some sequence parameters.
    params = pge2.check(psq, sys_ge, 'PNSwt', PNSwt);

    % Check accuracy of the psq sequence representation against the .seq file
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
    seq.read([seq_file_path seq_name '.seq']);
    xml_path = []; % if nonempty, compare against WTools/Pulse Studio output
    pge2.validate(psq, sys_ge, seq, xml_path, 'row', [], 'plot', false);

    % save to .mat file and add it to the tar archive
    pislquant = 2;  % only relevant for the 'adj_receive_gain.seq' sequence
    save([seq_name '.mat'], 'psq', 'params', 'pislquant');
    system(sprintf('tar --append -f %s %s.mat', tar_file_name, seq_name));

    fprintf('\n\n\n%s\n', repmat('-', 1, 79));

    % write to .pge file and add it to the tar archive
    pge2.serialize(psq, [seq_name '.pge'], 'pislquant', pislquant, 'params', params);
    system(sprintf('tar --append -f %s %s.pge', tar_file_name, seq_name));
end

% add .list file to tar archive
fclose(fid);
system(sprintf('tar --append -f %s %s', tar_file_name, scans_list_file));

