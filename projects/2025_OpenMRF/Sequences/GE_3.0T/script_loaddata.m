
% p = pge2.utils.normalizepath('/home/jon/dropbox/Shared/Pulseq/OpenMRF/20260212/Exam8521');
p = pge2.utils.normalizepath('/home/jon/transfer/Exam8521');

% pge722.entry, wasabi
S = [3 7 8 13];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
        save(['d_wasabi_series' num2str(s) '_run' num2str(ii)]);
        pause(2);
    end
end

% pge724.entry, mrf_yun
S = [4 9];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
        save(['d_mrf_yun_series' num2str(s) '_run' num2str(ii)]);
        pause(2);
    end
end

% pge726.entry, cmrf_t1_t2
S = [5 10];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
        save(['d_cmrf_t1_t2_series' num2str(s) '_run' num2str(ii)]);
        pause(2);
    end
end

% pge728.entry, cmrf_t1_t2_t1p_300Hz
S = [6 11];  % series numbers
for s = S
    pth = pge2.utils.normalizepath([p 'Series' num2str(s)]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for ii = 1:length(F)-1
        d = pge2.utils.loaddata([pth F{ii}]);
        save(['d_cmrf_t1_t2_g1p_300Hz_series' num2str(s) '_run' num2str(ii)]);
        pause(2);
    end
end

