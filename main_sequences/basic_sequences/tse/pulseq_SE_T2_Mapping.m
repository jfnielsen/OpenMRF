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
TE                  = [10; 15; 20; 30; 40; 60; 85; 120; 175; 250] *1e-3;
TSE_fixed.n_echo    = 1;
TSE_fixed.TR        = 5000 *1e-3;
TSE_fixed.exc_time  = 1.5  *1e-3;
TSE_fixed.rfc_time  = 2.5  *1e-3;
TSE_fixed.exc_tbw   = 4;
TSE_fixed.rfc_tbw   = 4;
TSE_fixed.t_acq     = 3.2 *1e-3 + 2*system.adcDeadTime;
TSE_fixed.os_mode   = 1;       % read oversampling: 0 off, 1 on
TSE_fixed.mode_exc  = 'sinc';  % 'sigpy_SLR' or 'sinc'
TSE_fixed.mode_rfc  = 'sinc';  % 'sigpy_SLR' or 'sinc'
TSE_fixed.enc_mode  = 'linear';

%% add sequence blocks
ndummy    = 1;
flag_plot = 1;

for loop_T2  = 1:numel(TE)
for loop_TR  = 1-ndummy : FOV.Ny
    TSE    = TSE_fixed;
    TSE.TE = TE(loop_T2);
    [TSE, ktraj_adc, ktraj_full] = TSE_init(TSE, FOV, system, flag_plot);
    TSE_add();
    ndummy    = 0;
    flag_plot = 0;
end
end

TSE.TE = TE;

%% plot sequence diagram
% seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();