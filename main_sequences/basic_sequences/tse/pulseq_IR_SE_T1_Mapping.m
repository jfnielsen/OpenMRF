%% init pulseq
% basic TSE (turbo-spin-echo) readout
clear
seq_name = 'IR_SE_T1_Mapping';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nx       = 128;
FOV.Ny       = 128;
FOV.fov_x    = 250 *1e-3;
FOV.fov_y    = 250 *1e-3;
FOV.dz       = 6 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% TSE sequence parameters
TSE.n_echo    = 1;
TSE.TE        = 7.3   *1e-3;
TSE.TR        = 5000  *1e-3;
TSE.exc_time  = 1.5 *1e-3;
TSE.rfc_time  = 2.5 *1e-3;
TSE.exc_tbw   = 4;
TSE.rfc_tbw   = 4;
TSE.t_acq     = 3.2 *1e-3 + 2*system.adcDeadTime;
TSE.os_mode   = 1;       % read oversampling: 0 off, 1 on
TSE.mode_exc  = 'sinc';  % 'sigpy_SLR' or 'sinc'
TSE.mode_rfc  = 'sinc';  % 'sigpy_SLR' or 'sinc'
TSE.enc_mode  = 'linear';

flag_plot = 1;
[TSE, ktraj_adc, ktraj_full] = TSE_init(TSE, FOV, system, flag_plot);

%% params: Inversion
INV.rf_type      = 'HYPSEC_inversion';
INV.tExc         = 10 *1e-3;  % [s]  hypsech pulse duration
INV.beta         = 700;       % [Hz] maximum rf peak amplitude
INV.mu           = 4.9;       % [ ]  determines amplitude of frequency sweep
INV.inv_rec_time = [35; 60; 95; 160; 260; 430; 710; 1180; 1940; 3200] *1e-3;
INV = INV_init(INV, FOV, system);

%% correct inversion recovery delays
dur_corr = mr.calcDuration([INV.gx_crush, INV.gy_crush, INV.gz_crush]) + mr.calcDuration(TSE.GS1) + mr.calcDuration(TSE.rf_exc)/2; % end of inversion <-> middle of 90° excitation
for j = 1:numel(INV.inv_rec_time)
    INV.inv_rec_delay(j,1) = mr.makeDelay(INV.inv_rec_time(j) - dur_corr);
end

%% add sequence blocks
ndummy = 1;
for loop_INV = 1:numel(INV.inv_rec_time)
for loop_TR  = 1-ndummy : TSE.nex
    INV_add();
    TSE_add();
    ndummy = 0;
end
end


%% plot sequence diagram
% seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();