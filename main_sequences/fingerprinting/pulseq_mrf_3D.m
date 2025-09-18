%% init pulseq
clear
seq_name = 'mrf_3D';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% select scanner
pulseq_scanner = 'Siemens_Vida_3T_MIITT';

% select pns sim orientation
pns_orientation = 'axial';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy      = 256;         % [ ] matrix size
FOV.Nz       = 48;          % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 96  *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% params: MRF flipangles and repetition times
MRF.pattern = 'sin_83';
MRF.FAs     = MRF_calc_FAs_sin([10, 30, 200; 1, 83, 200; 10, 10, 100]) *pi/180;
MRF.NR      = numel(MRF.FAs);
MRF.TRs     = 12.0 *1e-3 *ones(MRF.NR,1);

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
SPI.lim_gz_slew   = 0.65;        % [ ] reduce stimulation during slice excitation
SPI.lim_reph_slew = 0.65;        % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 2*FOV.Nz;   % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 2.8 *1e-3; % [s] time for spoiler and rewinder gradients
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

[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);
MRF.TRs = SPI.TR;

%% averages of center partitions
SPI.kz_indc  = find(SPI.kz_area==0);
SPI.kz_list  = [ 1:SPI.kz_indc-4, ...
                 SPI.kz_indc-3, ...
                 SPI.kz_indc-2, SPI.kz_indc-2, ...
                 SPI.kz_indc-1, SPI.kz_indc-1, SPI.kz_indc-1, ...
                 SPI.kz_indc,   SPI.kz_indc,   SPI.kz_indc,   SPI.kz_indc, ...
                 SPI.kz_indc+1, SPI.kz_indc+1, SPI.kz_indc+1, ...
                 SPI.kz_indc+2, SPI.kz_indc+2, ...
                 SPI.kz_indc+3, ...
                 SPI.kz_indc+4:FOV.Nz];

rng('default');
[~, temp_rand] = sort(rand(size(SPI.kz_list)));
SPI.kz_list = SPI.kz_list(temp_rand);

%% params: Saturation (magnetization reset)
SAT.mode            = 'on';
SAT.rf_type         = 'adiabatic_BIR4';
SAT.bir4_tau        = 10 *1e-3;  % [s]  bir4 pulse duration
SAT.bir4_f1         = 640;       % [Hz] maximum rf peak amplitude
SAT.bir4_beta       = 10;        % [ ]  am waveform parameter
SAT.bir4_kappa      = atan(10);  % [ ]  fm waveform parameter
SAT.bir4_dw0        = 30000;     % [rad/s] fm waveform scaling
SAT.sat_rec_time    = 3.0;       % [s] saturation recovery delay
SAT.crush_nTwists_z = FOV.Nz * 3.7;
SAT = SAT_init(SAT, FOV, system);

%% params: Inversion
INV.rf_type         = 'HYPSEC_inversion';
INV.tExc            = 10 *1e-3;  % [s]  hypsech pulse duration
INV.beta            = 700;       % [Hz] maximum rf peak amplitude
INV.mu              = 4.9;       % [ ]  determines amplitude of frequency sweep
INV.inv_rec_time    = [10] *1e-3;
INV.crush_nTwists_z = FOV.Nz * 2.9;
INV = INV_init(INV, FOV, system);

%% noise pre-scans
SPI.Nnoise = 64;
SPI_add_prescans();

%% create sequence

temp_start_stop = 1;

% pre-saturation
SAT_add();

for loop_kz = SPI.kz_indc%SPI.kz_list

    % start simulation of .seq file (kz=0)
    if (loop_kz == SPI.kz_indc) && (temp_start_stop==1)
        seq.addBlock(mr.makeLabel('SET', 'START', 1));
    end

    % saturation recovery
    SAT_add();

    % inversion
    INV_add();
    
    % spiral imaging
    for loop_NR = 1:SPI.NR
        SPI_add();
    end
    
    % stop simulation of .seq file (kz=0)
    if (loop_kz == SPI.kz_indc) && (temp_start_stop==1)
        seq.addBlock(mr.makeLabel('SET', 'STOP', 1));
        temp_start_stop = 0;
    end

end

%% plot sequence diagram
seq.plot('timeRange', [5 26]);

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();