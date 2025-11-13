%% init pulseq
clear
seq_name = 'IRFISP_MRF';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 1; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 1; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% select scanner
pulseq_scanner = 'United_Imaging_5T';

% select pns sim orientation
pns_orientation = 'coronal';

% init system, seq object and load pulseq user information
pulseq_init();

% maximum RF [Hz]
f1_max   = 800;
seq_name = [seq_name '_' num2str(f1_max) 'Hz'];

%% FOV geometry
FOV.Nxy      = 256;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 300  *1e-3;  % [m] FOV geometry
FOV.dz       = 5   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% params: MRF flipangles and repetition times
MRF.pattern         = 'yun'; % select pattern: e.g. 'yun' or 'cao'
[MRF.FAs, MRF.TRs]  = MRF_get_FAs_TRs(MRF.pattern, 1);
MRF.FAs(MRF.FAs==0) = 1e-6;
MRF.NR              = numel(MRF.FAs);

%% params: Spiral Readouts

% import params from MRF struct
SPI.NR             = MRF.NR;    % [ ] number of repetitions
SPI.mrf_import.TRs = MRF.TRs;   % [s] repetition times
SPI.mrf_import.FAs = MRF.FAs;   % [rad] flip angles

% slice excitation
SPI.exc_mode      = 'sinc';      % 'sinc' or 'sigpy_SLR'
SPI.exc_time      = 2.0 *1e-3;   % [s] excitation time
SPI.exc_tbw       = 6;           % [ ] time bandwidth product
SPI.exc_fa_mode   = 'import';    % 'equal',  'ramped',  'import'  
SPI.lim_gz_slew   = 0.5;         % [ ] reduce stimulation during slice excitation
SPI.lim_reph_slew = 0.5;         % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 4;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 2.0 *1e-3;  % [s] time for spoiler and rewinder gradients
SPI.lim_spoil_slew = 0.5;        % [ ] reduce stimulation during gradient spoiling

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
INV.tExc         = 8.0 *1e-3;  % [s]  hypsech pulse duration
INV.beta         = f1_max;     % [Hz] maximum rf peak amplitude
INV.mu           = 4.9;        % [ ]  determines amplitude of frequency sweep
INV.inv_rec_time = 0.01;       % [s]  inversion recovery time
INV = INV_init(INV, FOV, system);

%% noise pre-scans
SPI.Nnoise = 16;
SPI_add_prescans();

%% create sequence

% inversion
if flag_GE==1
	seq.addBlock(mr.makeDelay(system.gradRasterTime), mr.makeLabel('SET', 'TRID', 4));
end
INV_add();

% spiral imaging
for loop_NR = 1:SPI.NR
    SPI_add();
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();