%% init
% basis: spiral readout
% use for: B0 and B1 mapping via WASABI method -> doi.org/10.1002/mrm.26133
clear
seq_name = 'wasabi';

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

%% FOV geometry
FOV.Nxy      = 128;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 300  *1e-3;  % [m] FOV geometry
FOV.dz       = 8  *1e-3;    % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% params: Spiral Readouts

% basic params
SPI.NR     = 4;             % [ ] number of repetitions
SPI.TR     = 0 *1e-3;       % [s] repetition time, 0 for minimization
SPI.Ndummy = 0;             % [ ] initial dummy loops

% slice excitation
SPI.exc_mode       = 'sinc';        % 'sinc' or 'sigpy_SLR'
SPI.exc_shape      = 'ex';          % only for sigpy: 'st' or 'ex' 
SPI.exc_time       = 2.0 *1e-3;     % [s] excitation time
SPI.exc_tbw        = 4;             % [ ] time bandwidth product
SPI.exc_fa_mode    = 'equal';       % 'equal',  'ramped',  'import'  
SPI.exc_flipangle  = 90 *pi/180;    % [rad] const FA or start of FA ramp -> set [] for auto mode
SPI.lim_reph_slew  = 0.5;           % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 8;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 4.0 *1e-3;  % [s] time for spoiler and rewinder gradients
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
SPI.geo.Nvds     = 3.9;           % number of vds-spirals for sampling the kspce center
SPI.geo.BW       = 400 *1e3;      % [Hz] bandwidth of spiral acquisition
SPI.geo.Fcoeff   = [1  0];        % [1 0] for archimedean (equal density), [1 -0.99] for logarithmic (variable density)
SPI.geo.grad_lim = 1/sqrt(3);     % limit of gradient field strength
SPI.geo.slew_lim = 1/sqrt(3);     % limit of slew rate
SPI.geo.kmax     = SPI.kmax;      % determines resolution
SPI.geo.t_dwell  = 1/SPI.geo.BW;  % [s] dwell time for spiral acquisition

[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);
SPI.Trec = 1;

%% fat saturation and reset
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% WASABI preparation
WASABI.f0       = 127732435;     % [Hz]  larmor frequency for 3T
WASABI.ppm      = -2 : 0.1 : 2;  % [ppm] wasabi pulse offresonances
WASABI.f1       = 150;           % [Hz]  wasabi pulse amplitude
WASABI.tau      = 5 *1e-3;       % [s]   pulse duration (5ms)
WASABI.phase    = 0;             % [rad] pulse phase
WASABI.B0       = WASABI.f0 / system.gamma; % [T]
WASABI.B1       = WASABI.f1 / system.gamma; % [T]
WASABI.n_ppm    = numel(WASABI.ppm);        % [ ]
WASABI.f_off    = WASABI.ppm *1e-6 *WASABI.B0 *system.gamma; % [Hz] wasabi pulse offresonace frequencies
WASABI.tau      = round(WASABI.tau/system.rfRasterTime) * system.rfRasterTime; % prevent timing errors
WASABI          = WASABI_init(WASABI, FOV, system);

%% add sequence loops

ndummy = 1;
for loop_ppm = 1 : WASABI.n_ppm + 1
for loop_NR = 1-ndummy:SPI.NR

    ndummy = 0;

    % saturation, recovery, fat
    seq.addBlock(mr.makeLabel('SET', 'TRID', 4)); % TRID label    
    SAT_add(); % saturation    
    seq.addBlock(mr.makeDelay(SPI.Trec)); % recovery time    
    FAT_add(); % fat saturation

    % WASABI preparation
    if (loop_ppm < WASABI.n_ppm + 1)
        WASABI.rf.freqOffset = WASABI.f_off(loop_ppm);
        seq.addBlock(WASABI.rf);
    else
        WASABI.rf_dummy = WASABI.rf;
        WASABI.rf_dummy.signal = WASABI.rf_dummy.signal*0 + 1e-6;
        WASABI.rf_dummy.freqOffset = 0;
        seq.addBlock(WASABI.rf_dummy);
    end
    seq.addBlock(WASABI.gx_crush, WASABI.gy_crush, WASABI.gz_crush);

    % spiral readouts
    SPI_add();

end    
end

%% plot sequence diagram
seq.plot('TimeRange',[2 5]*SPI.Trec)

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();

%% export additional .seq file for receive gain adjustment
if flag_backup>0
    [seq_adj, external_path_adj] = GE_adj_receive_gain(system, 5, 2.0, SPI.adc, pi/2, FOV.dz, external_path, wip_id);
end