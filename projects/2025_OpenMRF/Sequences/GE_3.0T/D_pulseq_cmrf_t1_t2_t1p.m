%% init pulseq
clear
seq_name = 'cmrf_t1_t2_t1p';

% optional flags
flag_backup = 1; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 1; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 1; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% select scanner
% pulseq_scanner = 'Siemens_FreeMax_0,55T_MIITT';
% pulseq_scanner = 'Siemens_Sola_1,5T_MIITT';
% pulseq_scanner = 'Siemens_Vida_3T_MIITT';
pulseq_scanner = 'GE_Signa_3T_MIITT';

% select pns sim orientation
pns_orientation = 'coronal';

% init system, seq object and load pulseq user information
pulseq_init();

% store spin-lock amplitude to .seq name
fSL = 300;
seq_name = [seq_name '_' num2str(fSL) 'Hz'];

%% FOV geometry
FOV.Nxy      = 192;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 300  *1e-3;  % [m] FOV geometry
FOV.dz       = 8   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% params: MRF contrast encoding

% encoding list with contrast preparations for different segments
% 'No_Prep'    ->  use for recovery of longitudinal magnetization
% 'Saturation' ->  use for T1 encoding
% 'Inversion'  ->  use for T1 encoding
% 'T2'         ->  use for T2 encoding
% 'SL'         ->  use for T1p encoding
% 'ADIASL'     ->  use for adiabatic T1p encoding
% 'MLEV'       ->  use for T2 or T2p encoding

MRF.enc_list = {
'Inversion';
'No_Prep';
'T2';
'T2';
'Inversion';
'No_Prep';
'SL';
'SL';
'Inversion';
'No_Prep';
'T2';
'T2';
'Inversion';
'No_Prep';
'SL';
'SL'
};

MRF.n_segm = numel(MRF.enc_list);

%% params: MRF flipangles and repetition times
MRF.nr     = 48;                         % numer of readouts per hear beat
MRF.NR     = MRF.n_segm * MRF.nr;        % total number of readouts
MRF.TRs    = 0.0 *1e-3 *ones(MRF.NR,1);  % minimize TRs
MRF.FA_min = 4 *pi/180;                  % [rad] minimum flip angle
MRF.FA_max = 15 *pi/180;                 % [rad] minimum flip angle
MRF.FAs    = MRF_calc_FAs_sin_rand(MRF.FA_min, MRF.FA_max, MRF.nr, MRF.n_segm);

%% params: Spiral Readouts

% import params from MRF struct
SPI.nr             = MRF.nr;    % [ ] number of readouts per heart beat
SPI.NR             = MRF.NR;    % [ ] number of repetitions
SPI.mrf_import.TRs = MRF.TRs;   % [s] repetition times
SPI.mrf_import.FAs = MRF.FAs;   % [rad] flip angles

% slice excitation
SPI.exc_mode      = 'sinc';      % 'sinc' or 'sigpy_SLR'
SPI.exc_time      = 0.8 *1e-3;   % [s] excitation time
SPI.exc_tbw       = 2;           % [ ] time bandwidth product
SPI.exc_fa_mode   = 'import';    % 'equal',  'ramped',  'import' 
SPI.lim_gz_slew   = 0.9;         % [ ] reduce stimulation during slice excitation
SPI.lim_reph_slew = 0.9;         % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 4;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 0.8 *1e-3;  % [s] time for spoiler and rewinder gradients
SPI.lim_spoil_slew = 0.9;        % [ ] reduce stimulation during gradient spoiling

% rf spoiling
SPI.spoil_rf_mode = 'lin';      % rf spoiling mode: 'lin' or 'quad'
SPI.spoil_rf_inc  = 0 *pi/180;  % rf spoiling increment [rad]

% spiral geometry mode
SPI.geo.interleave_mode = 'RoundGoldenAngle';
SPI.geo.traj_mode       = 'vds';

% vds parameters
SPI.Nunique       = 48;            % number of unique projection angles
SPI.deltak        = 1/FOV.fov_xy;  % [1/m] kspace sampling
SPI.kmax          = SPI.deltak * FOV.Nxy/2;
SPI.geo.Nvds      = 24;            % number of vds-spirals for sampling the kspce center
SPI.geo.BW        = 500 *1e3;      % [Hz] bandwidth of spiral acquisition
SPI.geo.Fcoeff    = [1  -0.5];     % [1 0] for archimedean (equal density), [1 -0.5] for logarithmic (variable density)
SPI.geo.grad_lim  = 1/sqrt(3);     % limit of gradient field strength
SPI.geo.slew_lim  = 1/sqrt(3);     % limit of slew rate
SPI.geo.kmax      = SPI.kmax;      % determines resolution
SPI.geo.t_dwell   = 1/SPI.geo.BW;  % [s] dwell time for spiral acquisition

[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system, 1);
MRF.TRs = SPI.TR;

%% params: Inversion
INV.rf_type      = 'HYPSEC_inversion';
INV.tExc         = 10 *1e-3;  % [s]  hypsech pulse duration
INV.beta         = 700;       % [Hz] maximum rf peak amplitude
INV.mu           = 4.9;       % [ ]  determines amplitude of frequency sweep
INV.inv_rec_time = [0.01 35 380 130] *1e-3;
INV = INV_init(INV, FOV, system);

%% params: T2 preparation
T2.exc_mode   = 'adiabatic_BIR4';
T2.rfc_dur    = 2 *1e-3;   % [s]  duration of composite refocusing pulses
T2.bir4_tau   = 10 *1e-3;  % [s]  bir4 pulse duration
T2.bir4_f1    = 640;       % [Hz] maximum rf peak amplitude
T2.bir4_beta  = 10;        % [ ]  am waveform parameter
T2.bir4_kappa = atan(10);  % [ ]  fm waveform parameter
T2.bir4_dw0   = 30000;     % [rad/s] fm waveform scaling
T2.prep_times = [40 80 40 80] * 1e-3;  % [s] inversion times
T2            = T2_init(T2, FOV, system);

%% params: Spin-Lock

% spin-lock pulses
SL.relax_type = {'T1p'};                  % T1p or T2p or T2
SL.seq_type   = {'BSL'};                  % BSL or CSL or RESL
SL.tSL        = [40 80 40 80] *1e-3;      % [s]  SL time
SL.fSL        = fSL * ones(size(SL.tSL)); % [Hz] SL amplitude

% excitation pulses
SL.exc_mode  = 'adiabatic_AHP';  % 'adiabatic_AHP', 'sinc', 'sigpy_SLR' or 'bp'
SL.exc_time  = 3.0 *1e-3;        % [s] excitation time
SL.adia_wmax = 600 * 2*pi;       % [rad/s] amplitude of adiabatic pulse

% refocusing pulses
SL.rfc_mode = 'bp';              % 'bp', 'sinc', 'sigpy_SLR' or 'comp'
SL.rfc_time = 1.0 *1e-3;         % [s] refocusing time

SL = SL_init(SL, FOV, system);

%% params: Fat Saturation
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

%% check MRF encoding params
MRF_check_enc_list();

%% adjust dynamic segment delays
MRF_adjust_segment_delays();

%% Trigger Mode
% on:  Cardiac
% off: Abdominal or Phantom

MRF.mode_trig = 'off';

% calc fixed segment timings
MRF.acq_duration      = sum(SPI.TR(1:MRF.nr));
MRF.prep_acq_duration = MRF.prep_max + MRF.acq_duration;

if strcmp(MRF.mode_trig, 'on')
    MRF.delay_soft = mr.makeSoftDelay(0, 'acq_end', 'offset', -MRF.prep_acq_duration, 'factor', 1); % block_duration [s] = offset [s] + input [s] / factor
    TRIG_IN = mr.makeTrigger('physio1', 'system', system, 'delay', 10e-6, 'duration', 10e-3); % Input Trigger
else
    MRF.seg_duration = 1000 *1e-3; % [s] adjust segment duration
    MRF.delay_soft   = mr.makeDelay( round((MRF.seg_duration-MRF.prep_acq_duration)/system.gradRasterTime)*system.gradRasterTime ); % fixed delay
end

%% noise pre-scans
SPI.Nnoise = 16;
SPI_add_prescans();

%% create sequence
MRF_add_segments();

%% plot sequence diagram
seq.plot();

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();

%% export additional .seq file for receive gain adjustment
if flag_backup>0
    [seq_adj, external_path_adj] = GE_adj_receive_gain(system, 5, 2.0, SPI.adc, pi/2, FOV.dz, external_path, wip_id);
end