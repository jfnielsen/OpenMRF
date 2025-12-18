% Get Pulseq toolbox
system('git clone --branch v1.5.1 git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% Get toolbox to convert .seq file to a .pge file for execution on GE
system('git clone --branch main git@github.com:HarmonizedMRI/PulCeq.git');
addpath PulCeq/matlab
addpath PulCeq/matlab/DataHash

