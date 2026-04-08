
% p = pge2.utils.normalizepath('/home/jon/dropbox/Shared/Pulseq/OpenMRF/20260212/Exam8521');
p = normalizepath('/export/jfnielse/transfer/Exam8793');

% series numbers and correspond .entry numbers
S = [13 14 15 17 18];
opuser1 = [721 722 724 726 726];

% load raw data and save to .mat file
assert(numel(S) == numel(opuser1));
for ii = 1:length(S)
    pth = normalizepath([p 'Series' num2str(S(ii))]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    d = pge2.utils.loaddata([pth F{end-1}]);
    save(['d_' num2str(opuser1(ii))], '-v7.3');
    pause(2);
end
