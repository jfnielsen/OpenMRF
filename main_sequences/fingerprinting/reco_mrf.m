%% load mrf study
clear

study_path      = 'Q:/data/Pulseq/Rawdata/tomgr/Prisma/250912_OpenMRF_NIST/';
study_name_mrf  = 'meas_MID00571_FID167028_pulseq_mrf_rep1.dat';
study_name_traj = 'meas_MID00025_FID167054_pulseq_traj_mrf.dat';

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
                                [], ...                               % kz partitions for 3D MRF; used to eliminate unnecessary repetitions
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

%% optional: try reco
mixed_reco = mg_mixed_mrf_reco(DATA, ktraj(:,phi_id,:), Nxy, fov);

%% start MIITT reco

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

tic
% start reco & matching
[match, images] = mg_miitt_mrf_reco( DATA, ...              % [NCoils x NR x Nadc] mrf data
                                     ktraj, ...             % [2 x Nunique x Nadc] k-space trajectory
                                     fov, ...               % [m] field of view
                                     Nxy, ...               % [ ] matrix size
                                     phi_id, ...            % [NR x 1] projection phi ID
                                     DICT, ...              % [NR x Ndict] uncompressed dictionary OR structured field with compressed dictionary
                                     look_up, ...           % [Ndict x Nparams] look-up table
                                     look_up_names, ...     % {1 x Nparams} names of look-up table colums
                                     params_reco, ...       % parameters for basic recosntruction
                                     params_LR );           % parameters for low rank reconstruction
toc

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

%% optional: vis T1p match results
if isfield(match.LR, 'T1p')
    t1plims = t2lims;
    t1pcmp  = t2cmp;
    figure()
    ax1 = subplot(1,3,1);
        imagesc(match.direct.T1p, t1plims); axis image; axis off; colormap(gca, t1pcmp); colorbar;
        title('T1p direct');
    ax2 = subplot(1,3,2);
        imagesc(match.SVD.T1p, t1plims); axis image; axis off; colormap(gca, t1pcmp); colorbar;
        title('T1p SVD');
    ax3 = subplot(1,3,3);
        imagesc(match.LR.T1p, t1plims); axis image; axis off; colormap(gca, t1pcmp); colorbar;
        title('T1p LR');
    linkaxes([ax1 ax2 ax3]);
    clear ax1 ax2 ax3;
end

%% optional: vis db1 match results
if isfield(match.LR, 'db1')
    db1lims = [min(P.db1), max(P.db1)];
    db1cmp  = get_cmp('blue_red', 1000);
    figure()
    ax1 = subplot(1,3,1);
        imagesc(match.direct.db1, db1lims); axis image; axis off; colormap(gca, db1cmp); colorbar;
        title('dB1+ direct');
    ax2 = subplot(1,3,2);
        imagesc(match.SVD.db1, db1lims); axis image; axis off; colormap(gca, db1cmp); colorbar;
        title('dB1+ SVD');
    ax3 = subplot(1,3,3);
        imagesc(match.LR.db1, db1lims); axis image; axis off; colormap(gca, db1cmp); colorbar;
        title('dB1+ LR');
    linkaxes([ax1 ax2 ax3]);
    clear ax1 ax2 ax3;
end

%% optional: compare to NIST ref values
if 0
    [T1_ref, T2_ref] = NIST_references('1.5T_MnCl2');
    
    if ~exist('x', 'var')
        figure();
        imagesc(match.LR.T1, [0 3]); axis image; axis off; colormap(gca, t1cmp); colorbar;
        [x, y] = ginput(14);
    end
    
    r = 4; % radius of ROIs
    
    T1_mrf = zeros(14, 2);
    T2_mrf = zeros(14, 2);
    for j=1:14
        [tempx, tempy] = meshgrid(1:size(match.LR.M0,1), 1:size(match.LR.M0,2));
        roi = ((tempx - x(j)).^2 + (tempy - y(j)).^2) <= r^2;
        temp_t1 = match.LR.T1(roi);
        T1_mrf(j,1) = mean(temp_t1);
        T1_mrf(j,2) = std(temp_t1);
        temp_t2 = match.LR.T2(roi);
        T2_mrf(j,1) = mean(temp_t2);
        T2_mrf(j,2) = std(temp_t2);
        clear tempx tempy temp_t1 temp_t2;
    end
    
    figure()
    subplot(2,2,1)
    loglog(T1_ref, T1_mrf(:,1), '.', 'MarkerSize', 20)
    hold on
    temp_lim = [min([T1_ref; T1_mrf(:,1)])*0.9, max([T1_ref; T1_mrf(:,1)])*1.1];
    loglog(temp_lim, temp_lim, 'k--', 'LineWidth', 2)
    xlabel('T1 ref [s]')
    ylabel('T1 mrf [s]')
    xlim(temp_lim)
    ylim(temp_lim)
    set(gca, 'Fontsize', 12)
    axis square;
    
    subplot(2,2,2)
    loglog(T2_ref, T2_mrf(:,1), '.', 'MarkerSize', 20)
    hold on
    temp_lim = [min([T2_ref; T2_mrf(:,1)])*0.9, max([T2_ref; T2_mrf(:,1)])*1.1];
    loglog(temp_lim, temp_lim, 'k--', 'LineWidth', 2)
    xlabel('T2 ref [s]')
    ylabel('T2 mrf [s]')
    xlim(temp_lim)
    ylim(temp_lim)
    set(gca, 'Fontsize', 12)
    axis square;
    clear temp_lim;
    
    subplot(2,2,3)
    hold on
    plot(T1_ref, (T1_mrf(:,1)-T1_ref)./T1_ref*100, '.', 'MarkerSize', 20)
    yline(0, 'k--', 'LineWidth', 2)
    ylim([-20 20])
    xlabel('T1 ref')
    ylabel('T1 deviation [%]')
    set(gca, 'Fontsize', 12)
    axis square;
    
    subplot(2,2,4)
    hold on
    plot(T1_ref, (T2_mrf(:,1)-T2_ref)./T2_ref*100, '.', 'MarkerSize', 20)
    yline(0, 'k--', 'LineWidth', 2)
    ylim([-20 20])
    xlabel('T1 ref')
    ylabel('T2 deviation [%]')
    set(gca, 'Fontsize', 12)
    axis square;
end

%% save results
res.images = images;
res.match = match;
save_study_results(study_info, res);