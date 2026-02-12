
p = pge2.utils.normalizepath('/home/jon/dropbox/Shared/Pulseq/OpenMRF/20260212/Exam8521');

if 0

% pge722.entry, wasabi
S = [3 7 8 13];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
    end
    save(['d_wasabi_series' num2str(s) '_run' num2str(ii)]);
end

% pge724.entry, mrf_yun
S = [4 9];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
    end
    save(['d_mrf_yun_series' num2str(s) '_run' num2str(ii)]);
end
end

% pge726.entry, cmrf_t1_t2
S = [5 10];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
    end
    save(['d_cmrf_t1_t2_series' num2str(s) '_run' num2str(ii)]);
end

% pge728.entry, cmrf_t1_t2_t1p_300Hz
S = [6 11];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
    end
    save(['d_cmrf_t1_t2_g1p_300Hz_series' num2str(s) '_run' num2str(ii)]);
end

return

d = pge2.utils.loaddata([p 'ScanArchive_7347633TMRFIX_20260212_154428315.h5']);
save -v7.3 d722.mat d

return

d = pge2.utils.loaddata([p 'Series4/ScanArchive_7347633TMRFIX_20251218_182031072.h5']);
save -v7.3 d724.mat d

d = pge2.utils.loaddata([p 'Series5/ScanArchive_7347633TMRFIX_20251218_182136207.h5']);
save -v7.3 d726.mat d

d = pge2.utils.loaddata([p 'Series6/ScanArchive_7347633TMRFIX_20251218_182526493.h5']);
save -v7.3 d728.mat d

