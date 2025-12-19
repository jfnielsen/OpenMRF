% Get Pulseq toolbox
system('git clone --branch v1.5.1 git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab
warning('OFF', 'mr:restoreShape');  % turn off Pulseq warning for spirals

% Get toolbox to convert .seq file to a .pge file for execution on GE
system('git clone --branch main git@github.com:HarmonizedMRI/PulCeq.git');
addpath PulCeq/matlab
addpath PulCeq/matlab/DataHash

% GErecon function for loading ScanArchive files
addpath ~/Programs/orchestra-sdk-2.1-1.matlab/

