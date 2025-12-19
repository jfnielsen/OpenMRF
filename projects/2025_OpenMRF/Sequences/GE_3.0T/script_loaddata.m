
p = pge2.utils.normalizepath('/home/jon/transfer/20251218/Exam8318');

% pge722.entry, 251116_0307_mgram_wasabi
d = pge2.utils.loaddata([p 'Series3/ScanArchive_7347633TMRFIX_20251218_181648891.h5']);
save -v7.3 d722.mat d

% pge724.entry, 251116_0308_mgram_mrf_yun
d = pge2.utils.loaddata([p 'Series4/ScanArchive_7347633TMRFIX_20251218_182031072.h5']);
save -v7.3 d724.mat d

% pge726.entry, 251116_0309_mgram_cmrf_t1_t2
d = pge2.utils.loaddata([p 'Series5/ScanArchive_7347633TMRFIX_20251218_182136207.h5']);
save -v7.3 d726.mat d

% pge728.entry, 251116_0310_mgram_cmrf_t1_t2_t1p_300Hz
d = pge2.utils.loaddata([p 'Series6/ScanArchive_7347633TMRFIX_20251218_182526493.h5']);
save -v7.3 d728.mat d

