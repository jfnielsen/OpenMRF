
% p = pge2.utils.normalizepath('/home/jon/dropbox/Shared/Pulseq/OpenMRF/20260212/Exam8521');
p = normalizepath('/export/jfnielse/transfer/Exam8793');
p = normalizepath('/home/jon/transfer/Exam8797');

% series numbers and correspond .entry numbers
S = [11 12 13 14 15];
opuser1 = [721 722 724 725 726];

% load 5 most recent raw data files and save each to .mat file
assert(numel(S) == numel(opuser1));
for ii = 1:length(S)
    pth = normalizepath([p 'Series' num2str(S(ii))]);
    F = strsplit(ls(pth, '-tr'));   % .h5 file names
    for jj = 1 %:5
        d = pge2.utils.loaddata([pth F{end-jj}]);
        %fprintf('%s\n', [pth F{end-jj}]);
        save(sprintf('d_Series_%d_opuser1_%d_run_%d', S(ii), opuser1(ii), jj), 'd', '-v7.3');
    end
    pause(1);
end
