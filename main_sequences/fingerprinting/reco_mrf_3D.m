%% load mrf study
clear

study_path      = 'Q:/data/Pulseq/Rawdata/mgram/Skyra/250716_strukturphantom_git_test/';
study_name_mrf  = 'meas_MID00563_FID228481_pulseq_spi_flash_3d.dat';
% study_name_traj = 'meas_MID00035_FID92121_250725_1452_tomgr_traj_250725_1421_cor.dat';

% load twix_object, study info and pulseq meta data
[twix_obj, study_info, PULSEQ] = pulseq_read_meas_siemens([study_path study_name_mrf]);

% find path of the original .seq file
[~, ~, seq_path] = pulseq_get_user_definitions();
seq_path = [seq_path '/Pulseq_Workspace/' PULSEQ.pulseq_user '/' PULSEQ.seq_id(1:6) '/' PULSEQ.seq_id '/' PULSEQ.seq_name '.seq' ];

% select simulation mode and dictionary compression energy
sim_mode    = 'BLOCH'; % 'BLOCH' or 'EPG'
comp_energy = 0.9999;  % 0 for uncompressed dictionary

% enter soft delay input
softDelay = []; % only for cMRF -> twix_obj.hdr.Meas.adFree(7...) * 1e-3

%% define dictionary and look-up table

% T1, T2
P.T1.range = [0.01,  4]; P.T1.factor = 1.025;
P.T2.range = [0.001, 3]; P.T2.factor = 1.025;
P = MRF_get_param_dict(P, {'T2<T1'});
look_up       = [P.T1, P.T2];
look_up_names = {'T1', 'T2'};

% T1, T2, T1p
% P.T1.range  = [0.01,  4]; P.T1.factor  = 1.05;
% P.T2.range  = [0.001, 3]; P.T2.factor  = 1.05;
% P.T1p.range = [0.001, 3]; P.T1p.factor = 1.05;
% P = MRF_get_param_dict(P, {'T2<T1', 'T2<T1p', 'T1p<T1'});
% look_up       = [P.T1, P.T2, P.T1p];
% look_up_names = {'T1', 'T2', 'T1p'};

% T1, T2, B1+ correction
% P.T1.range  = [0.01,  4]; P.T1.factor = 1.025;
% P.T2.range  = [0.001, 3]; P.T2.factor = 1.025;
% P.db1.range = [0.8, 1.2]; P.db1.step  = 0.025;
% P = MRF_get_param_dict(P, {'T2<T1'});
% look_up       = [P.T1, P.T2, P.db1];
% look_up_names = {'T1', 'T2', 'db1'};

%% caclulate dictionary

% read .seq file
[SEQ, SIM] = MRF_read_seq_file( seq_path, ...                         % path of the original .seq file which was measured
                                twix_obj.hdr.Meas.lFrequency, ...     % f0; used for actual frequency offsets in fat suppression, CEST or WASABI modules
                                twix_obj.image.timestamp*0.0025, ...  % adc times stamps on a 2.5ms raster; used for correction of trigger delays
                                softDelay,...                         % soft delay input; used for correction of sequence timings
                                1, ...                                % kz partitions for 3D MRF; used to eliminate unnecessary repetitions
                                [], ...                               % echo mode; default: 'spiral_out'
                                1e-6, ...                             % raster time for the simulation 
                                0);                                   % flag_plot

% optional: vis ECG and sequence timings (works for 3T Skyra data)
% MRF_cardio_check_invivo_timings(SEQ, [study_path study_name_mrf], twix_obj, PULSEQ);

switch sim_mode
    case 'EPG'
        SIM  = MRF_sim_pre(SIM, P, [], 'EPG', 0, 0);
        DICT = MRF_sim_EPG(SIM, P, comp_energy);
    case 'BLOCH'
        NIso = 1000; % number of isochromats
        sfac = 2;    % factor for out-of-slice simulation
        z    = linspace(-1/2, 1/2, NIso)' *sfac *PULSEQ.FOV.dz;
        SIM  = MRF_sim_pre(SIM, P, z, 'BLOCH', 1, 0); 
        DICT = MRF_sim_BLOCH(SIM, P, z, [], comp_energy);
end

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
        error('wrong trajectory file!');
    end
    ktraj = ktraj_meas;
    clear ktraj_meas ktraj_hash1 ktraj_hash2;
end

%% load mrf rawdata;   optional: noise pre-whitening
if isfield(PULSEQ.SPI, 'Nnoise')
    Nnoise = PULSEQ.SPI.Nnoise;
else
    Nnoise = 0;
end
[DATA, ~, ~, ~, NOISE] = SPI_get_rawdata(twix_obj, Nnoise);
DATA = permute(DATA, [3,1,2]); % mrf rawdata
DATA(:,:,1:adcNPad)  = [];     % adc padding
ktraj(:,:,1:adcNPad) = [];     % adc padding

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
DATA = single(DATA);
DATA = permute(DATA, [2, 3, 1]);

% re-structure data
DATA = reshape(DATA, [Nseg, NR, NRead, NCoils]); % to do: check if reshape is correct!!!

% sort kz partitions and average repetitions
DATA_3D    = complex(zeros(Nz, NR, NRead, NCoils, 'single'));
kz_missing = [];
for j = 1:Nz
    ind = find(PULSEQ.SPI.kz_list==j);
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

%% start MIITT reco

% !!! replace this by a for loop which reconstructs all sub-slices !!!
% !!! store the results (match, images) in a 3d format !!!
nz = 10; % choose a sub-slice
DATA_sub = squeeze(DATA_3D_FFTz(nz, :,:,:));

% general reconstruction parameters
params_reco.DirectMatching     = false; % do reco via direct matching
params_reco.DirectMatching_SVD = true;  % do reco via SVD compression of the dictionary before matching
params_reco.LowRank            = true;  % do reco via iterative low-rank reconstruction
params_reco.ROVIR              = true;  % use ROVIR coils for outer FOV artifact suppression
params_reco.CoilComp           = true;  % use SVD Coil Compression
params_reco.ESPIRiT            = true;  % use ESPIRiT or openadapt for calculating cmaps
params_reco.rovir_thresh       = 2;     % automatic thresholding for ROVIR
params_reco.NCoils_v           = 8;     % virtual coils
params_reco.readOS             = 2;     % read oversampling factor
params_reco.NBlocks            = 12;    % number of blocks for pattern matching

% only for cMRF: control motion with sliding window reco
if isfield(PULSEQ.MRF, 'n_segm') && params_reco.DirectMatching
    params_reco.NSeg = PULSEQ.MRF.n_segm;
end

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

% start reco & matching
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

%% vis T1 T2 match results
t1lims = [0 2000] *1e-3;
t2lims = [0 1000]  *1e-3;
t1cmp  = get_cmp('T1', 1000, 1);
t2cmp  = get_cmp('T2', 1000, 1);

if isfield(match, 'direct')
    figure('Name','match results')
    ax1 = subplot(3,3,1);
        imagesc(abs(match.direct.M0)); axis image; axis off; colormap(gca, gray); colorbar;
        title('M0 direct');
    ax2 = subplot(3,3,2);
        imagesc(abs(match.SVD.M0)); axis image; axis off; colormap(gca, gray); colorbar;
        title('M0 SVD');
    ax3 = subplot(3,3,3);
        imagesc(abs(match.LR.M0)); axis image; axis off; colormap(gca, gray); colorbar;
        title('M0 LR');
    
    ax4 = subplot(3,3,4);
        imagesc(match.direct.T1, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
        title('T1 direct');
    ax5 = subplot(3,3,5);
        imagesc(match.SVD.T1, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
        title('T1 SVD');
    ax6 = subplot(3,3,6);
        imagesc(match.LR.T1, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
        title('T1 LR');
    
    ax7 = subplot(3,3,7);
        imagesc(match.direct.T2, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
        title('T2 direct');
    ax8 = subplot(3,3,8);
        imagesc(match.SVD.T2, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
        title('T2 SVD');
    ax9 = subplot(3,3,9);
        imagesc(match.LR.T2, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
        title('T2 LR');
    
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9]);
    clear ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9;
else
    figure('Name','match results')
    ax1 = subplot(3,2,1);
        imagesc(abs(match.SVD.M0)); axis image; axis off; colormap(gca, gray); colorbar;
        title('M0 SVD');
    ax2 = subplot(3,2,2);
        imagesc(abs(match.LR.M0)); axis image; axis off; colormap(gca, gray); colorbar;
        title('M0 LR');
    
    ax3 = subplot(3,2,3);
        imagesc(match.SVD.T1, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
        title('T1 SVD');
    ax4 = subplot(3,2,4);
        imagesc(match.LR.T1, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
        title('T1 LR');
    
    ax5 = subplot(3,2,5);
        imagesc(match.SVD.T2, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
        title('T2 SVD');
    ax6 = subplot(3,2,6);
        imagesc(match.LR.T2, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
        title('T2 LR');
    
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6]);
    clear ax1 ax2 ax3 ax4 ax5 ax6;
end