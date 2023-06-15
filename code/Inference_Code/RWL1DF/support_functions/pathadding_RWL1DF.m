% Add paths for RWL1-DF testing on images.
%
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated August 21, 2012. 
% 


comp_opt = 'Desk';

if strcmp(comp_opt, 'Desk')
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/Wavelab850/'))
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/YUV/'))
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/Noiselet_Code/'))
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/Wavelet_Code/'))
    addpath(genpath('~/Desktop/Versioned/2011_RWL1DCS/trunk/'))
    addpath(genpath('~/Desktop/Versioned/2011_RWL1DCS/trunk/MRIdata/'));
    rmpath /home/adam/Desktop/Versioned/CodeBaseRepo/trunk/Wavelab850/Books/WaveTour/WTCh06/
elseif strcmp(comp_opt, 'cortex')
    addpath(genpath('/homeold/adam/Versioned/CodeBaseRepo/trunk/Wavelab850/'))
    addpath(genpath('/homeold/adam/Versioned/CodeBaseRepo/trunk/YUV/'))
    addpath(genpath('/homeold/adam/Versioned/CodeBaseRepo/trunk/Noiselet_Code/'))
    addpath(genpath('/homeold/adam/Versioned/CodeBaseRepo/trunk/Wavelet_Code/'))
    addpath(genpath('/homeold/adam/Versioned/CodeBaseRepo/trunk/TFOCS_v1.0a/'))
    addpath(genpath('/homeold/adam/Versioned/2011_RWL1DCS/trunk/'))
    addpath(genpath('/homeold/adam/Versioned/2011_RWL1DCS/trunk/MRIdata/'));
    rmpath /homeold/adam/Versioned/CodeBaseRepo/trunk/Wavelab850/Books/WaveTour/WTCh06/
elseif strcmp(comp_opt, 'Lap')
    addpath(genpath('~/Versioned/CodeBaseRepo/trunk/Wavelab850/'))
    addpath(genpath('~/Versioned/2011_RWL1DCS/trunk/'))
    addpath(genpath('~/Versioned/CodeBaseRepo/trunk/SparsityToolboxes/toolbox_sparsity/toolbox/'))
    folder_name = '/home/adam/Versioned/2011_RWL1DCS/trunk/MRIdata/';
else
    error('Invalid computer!')
end

