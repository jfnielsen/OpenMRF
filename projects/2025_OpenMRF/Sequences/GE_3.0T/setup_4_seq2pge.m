% Get Pulseq toolbox
system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab
warning('OFF', 'mr:restoreShape');  % turn off Pulseq warning for spirals

% Get toolbox to convert .seq file to a PulSeg sequence (psq) object
system('git clone git@github.com:HarmonizedMRI/pulseg.git');
addpath pulseg/matlab
addpath(genpath('pulseg/matlab/third_party'));

% Get toolbox for plotting PulSeg (psq) object and exporting to binary file for GE
system('git clone git@github.com:HarmonizedMRI/pge2.git');
addpath pge2/matlab

% GErecon function for loading ScanArchive files
addpath ~/Programs/orchestra-sdk-2.1-1.matlab/

