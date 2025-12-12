% Create .pge files and put them in OpenMRF.tar

if ~exist('mr.Sequence');
    setup;
end

tarfn = 'OpenMRF.tar';  % put files here

seqdir = '~/dropbox/Pulseq/OpenMRF/';
seqdir = '~/Downloads/OpenMRF/';
seqdir = '~/Downloads/OpenMRF/GE_OpenMRF_251117_maybe_final/';
D = dir([seqdir '*.seq']);

CV1 = 721;   % for creating .entry file name

system(sprintf('rm -f %s', tarfn));
system(sprintf('tar cf %s main.m', tarfn));

for ii = 1:length(D)
    fn = replace(D(ii).name, '.seq', '');

    ceq = seq2ceq([seqdir fn '.seq'], 'usesRotationEvents', false);

    % scanner parameters/limits
    psd_rf_wait = 100e-6;  % RF-gradient delay, scanner specific (s)
    psd_grd_wait = 100e-6; % ADC-gradient delay, scanner specific (s)
    b1_max = 0.25;         % Gauss
    g_max = 5;             % Gauss/cm
    slew_max = 20;         % Gauss/cm/ms
    coil = 'xrm';          % MR750. See pge2.opts()
    sysGE = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

    % Check PNS and b1/gradient against scanner limits,
    % and extract some sequence parameters that we will pass to writeceq()
    pars = pge2.check(ceq, sysGE);

    % Plot the beginning of the sequence
    %S = pge2.plot(ceq, sysGE, 'timeRange', [0 0.02], 'rotate', false); 
    %pause(1);

    % Write ceq object to file.
    % pislquant is the number of ADC events used to set Rx gains in Auto Prescan
    writeceq(ceq, [fn '.pge'], 'pislquant', 2);

    % write entry file
    pge2.writeentryfile(CV1+ii-1, fn);

    % add files to tar file
    system(sprintf('tar --append -f %s pge%d.entry %s', tarfn, CV1+ii-1, [fn '.pge']));

    % After simulating in WTools/VM or scanning, grab the xml files 
    % and compare with the seq object:
    seq = mr.Sequence();
    seq.read([seqdir fn '.seq']);
    warning('OFF', 'mr:restoreShape');  % turn off Pulseq warning for spirals
    xmlPath = '~/transfer/xml/';
    pge2.validate(ceq, sysGE, seq, [], 'row', [], 'plot', false);
end
