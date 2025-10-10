%% load mrf study
clear

study_path       = '/Q/data/Pulseq/Rawdata/mgram/PrismaExpRad/251010_MRF_3D_Legs/';
study_name_mrf_1 = 'meas_MID00568_FID05672_mrf_3D_part1.dat';
study_name_mrf_2 = 'meas_MID00569_FID05673_mrf_3D_part2.dat';
study_name_mrf_3 = 'meas_MID00570_FID05674_mrf_3D_part3.dat';
study_name_mrf_4 = 'meas_MID00571_FID05675_mrf_3D_part4.dat';
study_name_traj  = 'meas_MID00743_FID05840_trajectory.dat';

% load twix_object, study info and pulseq meta data
[twix_obj_1, study_info_1, PULSEQ]   = pulseq_read_meas_siemens([study_path study_name_mrf_1]);
[twix_obj_2, study_info_2, PULSEQ_2] = pulseq_read_meas_siemens([study_path study_name_mrf_2]);
[twix_obj_3, study_info_3, PULSEQ_3] = pulseq_read_meas_siemens([study_path study_name_mrf_3]);
[twix_obj_4, study_info_4, PULSEQ_4] = pulseq_read_meas_siemens([study_path study_name_mrf_4]);

% load kz lists
kz_list_1 = PULSEQ.SPI.kz_list';
kz_list_2 = PULSEQ_2.SPI.kz_list';
kz_list_3 = PULSEQ_3.SPI.kz_list';
kz_list_4 = PULSEQ_4.SPI.kz_list';
clear PULSEQ_2 PULSEQ_3 PULSEQ_4;

% find path of the original .seq file
[~, ~, seq_path] = pulseq_get_user_definitions();
seq_path = [seq_path '/Pulseq_Workspace/' PULSEQ.pulseq_user '/' PULSEQ.seq_id(1:6) '/' PULSEQ.seq_id '/' PULSEQ.seq_name '.seq' ];

% select simulation mode and dictionary compression energy
sim_mode    = 'BLOCH'; % 'BLOCH' or 'EPG'
comp_energy = 0.9999;  % 0 for uncompressed dictionary

%% get sequence parameters
NR      = PULSEQ.SPI.NR;
Nxy     = PULSEQ.FOV.Nxy;
fov     = PULSEQ.FOV.fov_x;
ktraj   = PULSEQ.ktraj_reco;
phi_id  = PULSEQ.SPI.phi_id;
adcNPad = PULSEQ.SPI.adcNPad;

%% optional: use measured k-space trajectory
if exist('study_name_traj', 'var')
    ktraj_hash1               = pulseq_get_wave_hash(ktraj(:));
    [ktraj_meas, ktraj_hash2] = TRAJ_reco([study_path study_name_traj], 8:8:48);
    if ~strcmp(ktraj_hash1, ktraj_hash2)
        warning('wrong trajectory file!');
    end
    ktraj = ktraj_meas;
    clear ktraj_meas ktraj_hash1 ktraj_hash2;
end

%% load mrf rawdata; join parts; noise pre-whitening

% check number of noise pre-scans
if isfield(PULSEQ.SPI, 'Nnoise')
    Nnoise = PULSEQ.SPI.Nnoise;
else
    Nnoise = 0;
end

% load rawdata and noise pre-scans
[DATA_1, ~, ~, ~, NOISE] = SPI_get_rawdata(twix_obj_1, Nnoise); clear twix_obj_1;
DATA_2                   = SPI_get_rawdata(twix_obj_2, Nnoise); clear twix_obj_2;
DATA_3                   = SPI_get_rawdata(twix_obj_3, Nnoise); clear twix_obj_3;
DATA_4                   = SPI_get_rawdata(twix_obj_4, Nnoise); clear twix_obj_4;

% convert to single
DATA_1 = single(DATA_1);
DATA_2 = single(DATA_2);
DATA_3 = single(DATA_3);
DATA_4 = single(DATA_4);

% join DATA
DATA    = cat(1, DATA_1, DATA_2, DATA_3, DATA_4);
kz_list = [kz_list_1; kz_list_2; kz_list_3; kz_list_4];
clear DATA_1 DATA_2 DATA_3 DATA_4 kz_list_1 kz_list_2 kz_list_3 kz_list_4;

% permute mrf data
DATA = permute(DATA, [3,1,2]);

% adc padding
DATA(:,:,1:adcNPad)  = [];
ktraj(:,:,1:adcNPad) = [];

% noise pre-whitening
if ~isempty(NOISE)
    NOISE = permute(NOISE, [3,1,2]);  % permute noise prescans
    NOISE(:,:,1:adcNPad) = [];        % adc padding
    DATA = mg_noise_prewhitening(DATA, NOISE, 'cholesky', 1);
end
clear NOISE;

%% sort kz partitions and fft in z direction

% read dimensions, convert to single
[NCoils, Nseg, NRead] = size(DATA);
Nz   = PULSEQ.FOV.Nz;
Nseg = Nseg / NR;
DATA = permute(DATA, [2, 3, 1]);

% re-structure data
DATA = reshape(DATA, [NR, Nseg, NRead, NCoils]);
DATA = permute(DATA, [2, 1, 3, 4]);

% sort kz partitions and average repetitions
DATA_3D    = complex(zeros(Nz, NR, NRead, NCoils, 'single'));
kz_missing = [];
for j = 1:Nz
    ind = find(kz_list==j);
    if isempty(ind)
        warning(['missing k-space partition no.: ' num2str(j)]);
        kz_missing = [kz_missing; j]; % store missing partitions -> GRAPPA?
    else
        temp_data = DATA(ind, :,:,:);
        if size(temp_data,1) == 1
            temp_data = squeeze(temp_data);
        else
            temp_data = squeeze(mean(temp_data));
        end
        DATA_3D(j,:,:,:) = temp_data;
    end
end
clear j ind temp_data DATA;

% fft along kz direction
DATA_3D_FFTz = kspace2image(DATA_3D, [1, 0, 0, 0]);
DATA_3D_FFTz = permute(DATA_3D_FFTz, [1, 4, 2, 3]);
clear DATA_3D;

%% read .seq file
[SEQ, SIM] = MRF_read_seq_file( seq_path, ...     % path of the original .seq file which was measured
                                [],   ...         % f0; used for actual frequency offsets in fat suppression, CEST or WASABI modules
                                [],   ...         % adc times stamps on a 2.5ms raster; used for correction of trigger delays
                                [],   ...         % soft delay input; used for correction of sequence timings
                                1,    ...         % kz partitions for 3D MRF; used to eliminate unnecessary repetitions
                                [],   ...         % echo mode; default: 'spiral_out'
                                1e-6, ...         % raster time for the simulation 
                                0);               % flag_plot

%% define dictionary and look-up table
P.T1.range = [0.01,  4]; P.T1.factor = 1.025;
P.T2.range = [0.001, 3]; P.T2.factor = 1.025;
P = MRF_get_param_dict(P, {'T2<T1'});
look_up       = [P.T1, P.T2];
look_up_names = {'T1', 'T2'};

%% reconstruction and low rank parameters

% general reconstruction parameters
params_reco.DirectMatching     = false; % do reco via direct matching
params_reco.DirectMatching_SVD = false; % do reco via SVD compression of the dictionary before matching
params_reco.LowRank            = true;  % do reco via iterative low-rank reconstruction
params_reco.ROVIR              = true;  % use ROVIR coils for outer FOV artifact suppression
params_reco.CoilComp           = true;  % use SVD Coil Compression
params_reco.ESPIRiT            = true;  % use ESPIRiT or openadapt for calculating cmaps
params_reco.rovir_thresh       = 2;     % automatic thresholding for ROVIR
params_reco.NCoils_v           = 8;     % virtual coils
params_reco.readOS             = 2;     % read oversampling factor
params_reco.NBlocks            = 12;    % number of blocks for pattern matching

% Low Rank reconstruction parameters
params_LR = setupParameters_LowrankMRF2D( );
params_LR.numIter          = 50;          % max number of iterations    
params_LR.block_dim        = [6,6];       % locally low-rank patch size
params_LR.block_step       = 6;           % overlap between local patches (if equal to patch size above, then patches are non-overlapping)
params_LR.lambdaLLR        = 0.02;        % locally low-rank regularization
params_LR.lambdaSpatialTV  = 0.003;       % spatial TV regularization
params_LR.lambdaWav        = 0;           % wavelet regularization
params_LR.betaMethod       = 'Dai-Yuan';
params_LR.beta             = 0.6;
params_LR.alpha            = 0.01;
params_LR.stopThresh       = 0;
params_LR.updateFreq       = 0;

%% loop sub slices

M0_Maps = zeros(Nz, Nxy, Nxy);
T1_Maps = zeros(Nz, Nxy, Nxy);
T2_Maps = zeros(Nz, Nxy, Nxy);
IP_Maps = zeros(Nz, Nxy, Nxy);
PC1s    = zeros(Nz, Nxy, Nxy);

for nz = 1:Nz

DATA_sub = squeeze(DATA_3D_FFTz(nz, :,:,:));

% caclulate dictionary for sub-slice
switch sim_mode
    case 'EPG'
        SIM_temp = MRF_sim_pre(SIM, P, [], 'EPG', 0, 0);
        DICT     = MRF_sim_EPG(SIM_temp, P, comp_energy);
    case 'BLOCH'
        NIso     = 1000; % number of isochromats
        z        = PULSEQ.FOV.dz / PULSEQ.FOV.Nz * ( linspace(0, 1, NIso)' + nz - PULSEQ.FOV.Nz/2 - 1);
        SIM_temp = MRF_sim_pre(SIM, P, z, 'BLOCH', 1, 0); 
        DICT     = MRF_sim_BLOCH(SIM_temp, P, z, [], comp_energy);
end

% start MIITT reco
[match, images] = mg_miitt_mrf_reco( DATA_sub, ...          % [NCoils x NR x Nadc] mrf data
                                     ktraj, ...             % [2 x Nunique x Nadc] k-space trajectory
                                     fov, ...               % [m] field of view
                                     Nxy, ...               % [ ] matrix size
                                     phi_id, ...            % [NR x 1] projection phi ID
                                     DICT, ...              % [NR x Ndict] uncompressed dictionary OR structured field with compressed dictionary
                                     look_up, ...           % [Ndict x Nparams] look-up table
                                     look_up_names, ...     % {1 x Nparams} names of look-up table colums
                                     params_reco, ...       % parameters for basic recosntruction
                                     params_LR );           % parameters for low rank reconstruction

M0_Map = match.LR.M0;
T1_Map = match.LR.T1;
T2_Map = match.LR.T2;
IP_Map = match.LR.IP;
PC1    = squeeze(images.LR(1,:,:));
save(['/Q/home/mgram/results_3d_legs/subslice_' num2str(nz) '.mat'], 'M0_Map', 'T1_Map', 'T2_Map', 'IP_Map', 'PC1');

M0_Maps(nz,:,:) = M0_Map;
T1_Maps(nz,:,:) = T1_Map;
T2_Maps(nz,:,:) = T2_Map;
IP_Maps(nz,:,:) = IP_Map;
PC1s(nz,:,:)    = PC1;

end

save(['results_subslices.mat'], 'M0_Maps', 'T1_Maps', 'T2_Maps', 'IP_Maps', 'PC1s');

%% vis

% t1lims = [0 2000]/3 *1e-3;
% t2lims = [0 1000]/3  *1e-3;
% t1cmp  = get_cmp('T1', 1000, 1);
% t2cmp  = get_cmp('T2', 1000, 1);
% 
% mask = mg_get_mask_fit(abs(PC1), 'holes');
% 
% figure('Name','match results')
% ax2 = subplot(1,2,1);
% imagesc(T1_Map.*mask, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
% title('T1 LR');
% ax3 = subplot(1,2,2);
% imagesc(T2_Map.*mask, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
% title('T2 LR');