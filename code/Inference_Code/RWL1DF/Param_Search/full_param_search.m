%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calc_l1 = 1;
calc_rwl1 = 1;
calc_dl1 = 1;
calc_drwl1 = 1;
calc_dcsamp = 0;
calc_modcs = 1;
calc_rwl2 = 0;
calc_wl1p = 1;
calc_drwl1b = 1;

% Some simulation parameters
reset_samp = 0;
sim_type = 'com_chan';

sim_param_setup

% Generate Simulated CS measurements
[y,MEAS_FUN,x0] = Load_Data(sim_type,filename, load_params);

% General Parameters
TOL = 1e-3;          % TFOCS tolerance parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMSE_mins = zeros(1,7);

% BPDN search
cparamsearch_BPDN
params_BPDN = min_params;
RMSE_mean = mean(RMSE_KEEP, numel(size(RMSE_KEEP)));
RMSE_mins(1) = min(RMSE_mean(:));

% RWL1 search
cparamsearch_RWL1 
params_RWL1 = min_params;
RMSE_mean = mean(RMSE_KEEP, numel(size(RMSE_KEEP)));
RMSE_mins(2) = min(RMSE_mean(:));

% BPDN-DF search
cparamsearch_BPDNDF 
params_BPDNDF = min_params;
RMSE_mean = mean(RMSE_KEEP, numel(size(RMSE_KEEP)));
RMSE_mins(3) = min(RMSE_mean(:));

% RWL1-DF search
cparamsearch_RWL1DF 
params_RWL1DF = min_params;
RMSE_mean = mean(RMSE_KEEP, numel(size(RMSE_KEEP)));
RMSE_mins(4) = min(RMSE_mean(:));

% RWL1-DF search
cparamsearch_RWL1DFb
params_RWL1DFb = min_params;
RMSE_mean = mean(RMSE_KEEP, numel(size(RMSE_KEEP)));
RMSE_mins(5) = min(RMSE_mean(:));

% RWL1-DF search
cparamsearch_WL1DF
params_WL1DF = min_params;
RMSE_mean = mean(RMSE_KEEP, numel(size(RMSE_KEEP)));
RMSE_mins(6) = min(RMSE_mean(:));

% RWL1-DF search
cparamsearch_MODCS
params_MODCS = min_params;
RMSE_mean = mean(RMSE_KEEP, numel(size(RMSE_KEEP)));
RMSE_mins(7) = min(RMSE_mean(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
