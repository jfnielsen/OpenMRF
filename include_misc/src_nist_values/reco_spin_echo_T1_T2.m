% reconstruction of Siemens T1 and T2 weighted spin-echo images
% use for T1/T2 mapping of the NIST phantom
% mgram; V1; 14.09.2025

%% start reco
clear

study = 4;

switch study

    case 1  % rawdata: Avanto 1.5T (Wuerzburg)    
    ref_name       =  '1.5T_MnCl2_Wuerzburg';
    study_path     = 'Q:/data/Pulseq/Rawdata/tomgr/Avanto/250912_OpenMRF_NIST';
    study_names_t1 = { 'meas_MID00845_FID16202_se_ir_35.dat',
                       'meas_MID00846_FID16203_se_ir_60.dat',
                       'meas_MID00847_FID16204_se_ir_95.dat',
                       'meas_MID00848_FID16205_se_ir_160.dat',
                       'meas_MID00849_FID16206_se_ir_260.dat',
                       'meas_MID00850_FID16207_se_ir_430.dat',
                       'meas_MID00851_FID16208_se_ir_710.dat',
                       'meas_MID00852_FID16209_se_ir_1180.dat',
                       'meas_MID00853_FID16210_se_ir_1940.dat',
                       'meas_MID00854_FID16211_se_ir_3200.dat'};
    study_names_t2 = { 'meas_MID00855_FID16212_se_te_10.dat',
                       'meas_MID00856_FID16213_se_te_15.dat',
                       'meas_MID00857_FID16214_se_te_20.dat',
                       'meas_MID00858_FID16215_se_te_30.dat',
                       'meas_MID00859_FID16216_se_te_40.dat',
                       'meas_MID00860_FID16217_se_te_60.dat',
                       'meas_MID00861_FID16218_se_te_85.dat',
                       'meas_MID00862_FID16219_se_te_120.dat',
                       'meas_MID00863_FID16220_se_te_175.dat',
                       'meas_MID00864_FID16221_se_te_250.dat'};

    case 2  % rawdata: Prisma 3.0T (Wuerzburg)    
    ref_name       =  '3.0T_MnCl2_Wuerzburg';
    study_path     = 'Q:/data/Pulseq/Rawdata/tomgr/Prisma/250912_OpenMRF_NIST';
    study_names_t1 = { 'meas_MID00544_FID167001_se_ir_35ms.dat',
                       'meas_MID00545_FID167002_se_ir_60ms.dat',
                       'meas_MID00546_FID167003_se_ir_95ms.dat',
                       'meas_MID00547_FID167004_se_ir_160ms.dat',
                       'meas_MID00548_FID167005_se_ir_260ms.dat',
                       'meas_MID00549_FID167006_se_ir_430ms.dat',
                       'meas_MID00550_FID167007_se_ir_710ms.dat',
                       'meas_MID00551_FID167008_se_ir_1180ms.dat',
                       'meas_MID00552_FID167009_se_ir_1940ms.dat',
                       'meas_MID00553_FID167010_se_ir_3200ms.dat'};
    study_names_t2 = { 'meas_MID00554_FID167011_se_te_10ms.dat',
                       'meas_MID00555_FID167012_se_te_15ms.dat',
                       'meas_MID00556_FID167013_se_te_20ms.dat',
                       'meas_MID00557_FID167014_se_te_30ms.dat',
                       'meas_MID00558_FID167015_se_te_40ms.dat',
                       'meas_MID00559_FID167016_se_te_60ms.dat',
                       'meas_MID00560_FID167017_se_te_85ms.dat',
                       'meas_MID00561_FID167018_se_te_120ms.dat',
                       'meas_MID00562_FID167019_se_te_175ms.dat',
                       'meas_MID00563_FID167020_se_te_250ms.dat'};

    case 3  % rawdata: Aera 1.5T (Freiburg)    
    ref_name       =  '1.5T_MnCl2_Freiburg';
    study_path     = 'Q:/data/Pulseq/Rawdata/tomgr/Aera/250906_OpenMRF_NIST';
    study_names_t1 = { 'meas_MID00076_FID24830_t1_se_cor_TI_35ms.dat',
                       'meas_MID00077_FID24831_t1_se_cor_TI_60ms.dat',
                       'meas_MID00078_FID24832_t1_se_cor_TI_95ms.dat',
                       'meas_MID00079_FID24833_t1_se_cor_TI_160ms.dat',
                       'meas_MID00080_FID24834_t1_se_cor_TI_260ms.dat',
                       'meas_MID00081_FID24835_t1_se_cor_TI_430ms.dat',
                       'meas_MID00082_FID24836_t1_se_cor_TI_710ms.dat',
                       'meas_MID00083_FID24837_t1_se_cor_TI_1180ms.dat',
                       'meas_MID00084_FID24838_t1_se_cor_TI_1940ms.dat',
                       'meas_MID00085_FID24839_t1_se_cor_TI_3200ms.dat'};
    study_names_t2 = { 'meas_MID00086_FID24840_t2_se_cor_TE_10ms.dat',
                       'meas_MID00087_FID24841_t2_se_cor_TE_15ms.dat',
                       'meas_MID00088_FID24842_t2_se_cor_TE_20ms.dat',
                       'meas_MID00089_FID24843_t2_se_cor_TE_30ms.dat',
                       'meas_MID00090_FID24844_t2_se_cor_TE_40ms.dat',
                       'meas_MID00091_FID24845_t2_se_cor_TE_60ms.dat',
                       'meas_MID00092_FID24846_t2_se_cor_TE_85ms.dat',
                       'meas_MID00093_FID24847_t2_se_cor_TE_120ms.dat',
                       'meas_MID00094_FID24848_t2_se_cor_TE_175ms.dat',
                       'meas_MID00095_FID24849_t2_se_cor_TE_250ms.dat'};

    case 4  % rawdata: CimaX 3.0T (Freiburg)    
    ref_name       =  '3.0T_MnCl2_Freiburg';
    study_path     = 'Q:/data/Pulseq/Rawdata/tomgr/CimaX/250823_OpenMRF_NIST';
    study_names_t1 = { 'meas_MID00373_FID31834_t1_se_cor_TI_35ms.dat',
                       'meas_MID00374_FID31835_t1_se_cor_TI_60ms.dat',
                       'meas_MID00375_FID31836_t1_se_cor_TI_95ms.dat',
                       'meas_MID00376_FID31837_t1_se_cor_TI_160ms.dat',
                       'meas_MID00377_FID31838_t1_se_cor_TI_260ms.dat',
                       'meas_MID00378_FID31839_t1_se_cor_TI_430ms.dat',
                       'meas_MID00379_FID31840_t1_se_cor_TI_710ms.dat',
                       'meas_MID00380_FID31841_t1_se_cor_TI_1180ms.dat',
                       'meas_MID00381_FID31842_t1_se_cor_TI_1940ms.dat',
                       'meas_MID00382_FID31843_t1_se_cor_TI_3200ms.dat'};
    study_names_t2 = { 'meas_MID00383_FID31844_t2_se_cor_TE_10ms.dat',
                       'meas_MID00384_FID31845_t2_se_cor_TE_15ms.dat',
                       'meas_MID00385_FID31846_t2_se_cor_TE_20ms.dat',
                       'meas_MID00386_FID31847_t2_se_cor_TE_30ms.dat',
                       'meas_MID00387_FID31848_t2_se_cor_TE_40ms.dat',
                       'meas_MID00388_FID31849_t2_se_cor_TE_60ms.dat',
                       'meas_MID00389_FID31850_t2_se_cor_TE_85ms.dat',
                       'meas_MID00390_FID31851_t2_se_cor_TE_120ms.dat',
                       'meas_MID00391_FID31852_t2_se_cor_TE_175ms.dat',
                       'meas_MID00392_FID31853_t2_se_cor_TE_250ms.dat'};

end

%% reconstruct images

% load rawdata from twix objects
for j=1:numel(study_names_t1)
    twix_obj = mapVBVD(fullfile(study_path, study_names_t1{j}), 'ignoreSeg', 'removeOS');
    twix_obj = twix_obj{end};
    twix_obj = squeeze(twix_obj.image());
    rawdata_T1(j,:,:,:) = permute(twix_obj,[2,3,1]);
    clear twix_obj;
end
for j=1:numel(study_names_t2)
    twix_obj = mapVBVD(fullfile(study_path, study_names_t2{j}), 'ignoreSeg', 'removeOS');
    twix_obj = twix_obj{end};
    twix_obj = squeeze(twix_obj.image());
    rawdata_T2(j,:,:,:) = permute(twix_obj,[2,3,1]);
    clear twix_obj;
end

% ifft
Images_coils_T1 = kspace2image(rawdata_T1);
Images_coils_T2 = kspace2image(rawdata_T2);
clear rawdata_T1 rawdata_T2;

% calculate cmaps via espirit
[cmaps_T1, ~, ~] = mg_espirit_cmaps(squeeze(Images_coils_T1(end,:,:,:)), [], [], [], []);  % use the first image as the reference
[cmaps_T2, ~, ~] = mg_espirit_cmaps(squeeze(Images_coils_T2(1,:,:,:)),   [], [], [], []);  % use the last image as the reference

% calculate coil combined images
for j=1:size(Images_coils_T1, 1)
    Images_T1(j,:,:,:) = squeeze(sum(squeeze(Images_coils_T1(j,:,:,:)) .* conj(cmaps_T1)));
end
for j=1:size(Images_coils_T2, 1)
    Images_T2(j,:,:,:) = squeeze(sum(squeeze(Images_coils_T2(j,:,:,:)) .* conj(cmaps_T2)));
end
clear j Images_coils_T1 Images_coils_T2 cmaps_T1 cmaps_T2;

% zero interpolation filling
zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;
Images_T1 = mg_zero_filling(Images_T1, zero_params);
Images_T2 = mg_zero_filling(Images_T2, zero_params);

% rotate images to real axis
Images_T1 = real(Images_T1 .* exp(-1i*repmat(angle(Images_T1(end,:,:)),[size(Images_T1,1),1,1]))); % use the last image as the reference
Images_T2 = real(Images_T2 .* exp(-1i*repmat(angle(Images_T2(1,:,:)),[size(Images_T2,1),1,1])));   % use the first image as the reference

%% try to load mask_fit and rois
if isfile([ref_name '_mask_rois.mat'])
    load([ref_name '_mask_rois.mat']);
else
    mask_fit = mg_get_mask_fit(squeeze(mean(abs(Images_T1)))+squeeze(mean(abs(Images_T2))), 'holes');
end

%% T1 and T2 mapping
TI = [35 60 95 160 260 430 710 1180 1940 3200] *1e-3;
[T1_Map, ~, ~, R2_Map_T1] = mg_map_T1( double(real(Images_T1)), TI, mask_fit);

TE = [10 15 20 30 40 60 85 120 175 250] *1e-3;
[T2_Map, ~, R2_Map_T2] = mg_map_T12p( double(real(Images_T2)), TE, mask_fit );

%% select roi centers
t1cmp  = get_cmp('T1', 1000, 1);
t2cmp  = get_cmp('T2', 1000, 1);

[xlims, ylims] = meshgrid(1:size(mask_fit,1), 1:size(mask_fit,2));
xlims          = xlims(mask_fit==1);
ylims          = ylims(mask_fit==1);
xlims          = [min(xlims), max(xlims)] + 3*[-1 1];
ylims          = [min(ylims), max(ylims)] + 3*[-1 1];

if ~exist('xc', 'var')
    figure;
    ax1 = subplot(1,2,1);
    imagesc(T1_Map, [0 0.4]); axis image; axis off; colormap(ax1, t1cmp); xlim(xlims); ylim(ylims);
    ax2 = subplot(1,2,2);
    imagesc(T1_Map, [0 3]); axis image; axis off; colormap(ax2, t1cmp); xlim(xlims); ylim(ylims);
    [xc, yc] = ginput(13);
    save([ref_name '_mask_rois.mat'], 'xc', 'yc', 'mask_fit');
end

%% calc mean in ROIs
r   = 3;
roi = zeros(size(T2_Map));
for p=1:14
    [tempx, tempy] = meshgrid(1:size(T2_Map,1), 1:size(T2_Map,2));
    temp_roi = ((tempx - xc(p)).^2 + (tempy - yc(p)).^2) <= r^2;
    roi = roi + temp_roi*p;    
end
clear tempx tempy temp_roi;

T1 = zeros(14, 1);
T2 = zeros(14, 1);
for p=1:14
    temp_vals = T1_Map(roi==p);
    temp_r2   = R2_Map_T1(roi==p);
    temp_vals = temp_vals(temp_r2>0.975);
    T1(p,1)   = mean(temp_vals);
    temp_vals = T2_Map(roi==p);
    temp_r2   = R2_Map_T2(roi==p);
    temp_vals = temp_vals(temp_r2>0.95);
    T2(p,1)   = mean(temp_vals);
end
clear p temp_vals temp_r2;

%% vis results
figure;
ax1 = subplot(2,2,1);
imagesc(T1_Map, [0 0.4]); axis image; axis off; colormap(ax1, t1cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);
ax2 = subplot(2,2,2);
imagesc(T1_Map, [0 3]); axis image; axis off; colormap(ax2, t1cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);
ax3 = subplot(2,2,3);
imagesc(T2_Map, [0 0.05]); axis image; axis off; colormap(ax3, t2cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);
ax4 = subplot(2,2,4);
imagesc(T2_Map, [0 2]); axis image; axis off; colormap(ax4, t2cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);

%% compare to nominal values
[T1_nom, T2_nom] = NIST_references(ref_name(1:10));

figure()
subplot(1,2,1)
loglog(T1_nom, T2_nom, 'go', 'MarkerSize', 10)
hold on
loglog(T1, T2, 'k.', 'MarkerSize', 20)
xlabel('T1 [s]')
ylabel('T2 [s]')
legend('nom', 'meas', 'Location', 'best')
axis square
grid on
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'LineWidth', 2, 'Fontsize', 16); 

subplot(1,2,2)
semilogx(T1_nom, (T1-T1_nom)./T1_nom*100, 'r.', 'MarkerSize', 20)
hold on
semilogx(T2_nom, (T2-T2_nom)./T2_nom*100, 'b.', 'MarkerSize', 20)
yline(0, 'k--')
xlabel('T1/T2 [s]')
ylabel('deviation [%]')
legend('T1', 'T2', 'Location', 'best')
ylim([-50 50])
axis square
grid on
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'LineWidth', 2, 'Fontsize', 16); 

%% save results
save([ref_name '.mat'], 'T1_Map', 'T2_Map', 'R2_Map_T1', 'R2_Map_T2', 'T1', 'T2', 'roi');