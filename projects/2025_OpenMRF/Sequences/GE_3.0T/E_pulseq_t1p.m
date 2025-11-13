%% init pulseq
% basis: SPI readout
% use for: T1rho, T2rho or T2 mapping
clear
seq_name = 't1p_mapping';

% optional flags
flag_backup = 1; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 1; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% select scanner
% pulseq_scanner = 'Siemens_FreeMax_0,55T_MIITT';
% pulseq_scanner = 'Siemens_Sola_1,5T_MIITT';
% pulseq_scanner = 'Siemens_Vida_3T_MIITT';
pulseq_scanner = 'GE_Signa_3T_MIITT';

% select pns sim orientation
% pns_orientation = 'coronal';

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

%% params: Spiral Readouts

% basic params
SPI.NR     = 8;             % [ ] number of repetitions
SPI.TR     = 0 *1e-3;       % [s] repetition time, 0 for minimization
SPI.Ndummy = 0;             % [ ] initial dummy loops

% slice excitation
SPI.exc_mode       = 'sinc';        % 'sinc' or 'sigpy_SLR'
SPI.exc_time       = 2.0 *1e-3;     % [s] excitation time
SPI.exc_tbw        = 4;             % [ ] time bandwidth product
SPI.exc_fa_mode    = 'equal';       % 'equal',  'ramped',  'import'  
SPI.exc_flipangle  = 90 *pi/180;    % [rad] const FA or start of FA ramp -> set [] for auto mode
SPI.lim_gz_slew    = 0.8;           % [ ] reduce stimulation during slice excitation
SPI.lim_reph_slew  = 0.6;           % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 8;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 3.5 *1e-3;  % [s] time for spoiler and rewinder gradients
SPI.lim_spoil_slew = 0.5;        % [ ] reduce stimulation during gradient spoiling

% rf spoiling
SPI.spoil_rf_mode = 'quad';       % rf spoiling mode: 'lin' or 'quad'
SPI.spoil_rf_inc  = 117 *pi/180;  % rf spoiling increment [rad]

% spiral geoemtry mode
SPI.geo.traj_mode       = 'vds';         % 'standard' or 'vds'
SPI.geo.interleave_mode = 'Equal2Pi';    % 'Equal2Pi', 'Random2Pi', 'GoldenAngle', 'RandomGoldenAngle', 'RoundGoldenAngle'
SPI.deltak              = 1/FOV.fov_xy;  % [1/m] kspace sampling
SPI.kmax                = SPI.deltak * FOV.Nxy/2;

% spiral geometry: vds-hargreaves-toolbox
SPI.geo.Nvds     = 7.9;           % number of vds-spirals for sampling the kspce center
SPI.geo.BW       = 400 *1e3;      % [Hz] bandwidth of spiral acquisition
SPI.geo.Fcoeff   = [1  0];        % [1 0] for archimedean (equal density), [1 -0.99] for logarithmic (variable density)
SPI.geo.grad_lim = 1/sqrt(3);     % limit of gradient field strength
SPI.geo.slew_lim = 1/sqrt(3);     % limit of slew rate
SPI.geo.kmax     = SPI.kmax;      % determines resolution
SPI.geo.t_dwell  = 1/SPI.geo.BW;  % [s] dwell time for spiral acquisition

[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);

SPI.Trec = 5;

%% params: Spin-Lock

% spin-lock pulses
SL.relax_type = {'T1p'};                  % T1p or T2p or T2
SL.seq_type   = {'BSL'};                  % BSL or CSL or RESL
SL.tSL        = [8 : 8 : 80] *1e-3;       % [s]  SL time
SL.fSL        = fSL * ones(size(SL.tSL)); % [Hz] SL amplitude

% excitation pulses
SL.exc_mode  = 'adiabatic_AHP';  % 'adiabatic_AHP', 'sinc', 'sigpy_SLR' or 'bp'
SL.exc_time  = 3.0 *1e-3;        % [s] excitation time
SL.adia_wmax = 600 * 2*pi;       % [rad/s] amplitude of adiabatic pulse

% refocusing pulses
SL.rfc_mode = 'bp';              % 'bp', 'sinc', 'sigpy_SLR' or 'comp'
SL.rfc_time = 1.0 *1e-3;         % [s] refocusing time

SL = SL_init(SL, FOV, system);

% store spin-lock amplitude to .seq name
fSL = 300;
seq_name = [seq_name '_' num2str(fSL) 'Hz'];

%% fat saturation and reset
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_SL = 1 : SL.nSL    
for loop_NR = 1-SPI.Ndummy : SPI.NR

    % saturtion, recovery, fat
    seq.addBlock(mr.makeLabel('SET', 'TRID', 4)); % TRID label
    SAT_add(); % saturation
    seq.addBlock(mr.makeDelay(SPI.Trec)); % recovery time
    FAT_add(); % fat saturation

    % sl preparation
    seq.addBlock(mr.makeLabel('SET', 'TRID', 4+loop_SL)); % TRID label
    SL_add();

    % spiral imaging
    SPI_add();  

end    
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();

%% export additional .seq file for receive gain adjustment
if flag_backup>0
    [seq_adj, external_path_adj] = GE_adj_receive_gain(system, 5, 2.0, SPI.adc, pi/2, FOV.dz, external_path, wip_id);
end