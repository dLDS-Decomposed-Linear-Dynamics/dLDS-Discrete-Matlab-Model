%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Script that tests dictionary learning for BPDNDF dynamics %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath(genpath('/home/adam/GITrepos/pop_spike_dyn/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set parameters
block_size = 8;

sig_opts.N       = block_size^2;
sig_opts.M       = block_size^2;
sig_opts.S       = 2;
sig_opts.sig_var = 0.05;
sig_opts.sig_var = 1;
sig_opts.T_s     = 20;
sig_opts.nF      = 5;
sig_opts.sF      = 2;
sig_opts.dsamp   = 4;

sig_opts.noise_opts.p      = 0;
sig_opts.noise_opts.dg_var = 1e-4;

inf_opts.lambda_val      = 0.1;
inf_opts.lambda_history  = 0.90;
inf_opts.lambda_b        = 0.01;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.max_iter2       = 300;
inf_opts.special = '';

[D, F] = bpdndf_dynamics_learning(sig_opts, [], [], inf_opts);             % Run the dictionary learning algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a new dataset

block_size = 8;

sig_opts.N        = block_size^2;
sig_opts.M        = block_size^2;
sig_opts.S        = 2;
sig_opts.sig_var  = 0.05;
sig_opts.sig_var  = 1;
sig_opts.T_s      = 50000;
sig_opts.nF       = 2;
sig_opts.sF       = 2;
sig_opts.dsamp    = 4;
sig_params.unit_b = true;

sig_opts.noise_opts.p      = 0;
sig_opts.noise_opts.dg_var = 1e-4;

F_true = rand_dyn_create(sig_opts.M, sig_opts.nF, 'perm');
D_true = eye(sig_opts.M);

[X_ex, X_true, ~] = rand_seq_create(sig_opts, sig_opts.noise_opts, F_true);% Synthesize a new set of data examples
X_ex = D_true*X_ex; 
X_ex = X_ex(:,300:end);

%% run LDS
 
% Xdata       = iddata(X_ex.',zeros(size(X_ex,2),1),1);
% sys_builtin = ssest(Xdata,64);

%% Least Squares
imagesc(X_ex(:,2:end)/X_ex(:,1:(end-1)))


%% Dictionary Learning

block_size = 8;

inf_opts.N               = block_size^2;
% inf_opts.M               = block_size^2;
inf_opts.lambda_val      = 0.01;
inf_opts.lambda_history  = 0.1;
inf_opts.lambda_b        = 0.001;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.N_ex            = 150;
inf_opts.nF              = 2;
inf_opts.max_iter        = 10e3;
inf_opts.max_iter2       = 300;
inf_opts.step_d          = 50;
inf_opts.step_f          = 20;
inf_opts.plot_option     = 20;

inf_opts.special         = 'switched';
[D, F] = bpdndf_dynamics_learning(10*X_ex, [], [], inf_opts);              % Run the dictionary learning algorithm

%% BBC standard sysID

block_size = 12;

inf_opts.N               = 4*block_size^2;
% inf_opts.M               = block_size^2;
inf_opts.lambda_val      = 0.8;
inf_opts.lambda_history  = 0.9;
inf_opts.lambda_b        = 0.02;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.N_ex            = 150;
inf_opts.nF              = 1;
inf_opts.max_iter        = 10e3;
inf_opts.max_iter2       = 300;
inf_opts.step_d          = 50;
inf_opts.step_f          = 20;
inf_opts.plot_option     = 20;

inf_opts.special         = 'sysop';
[D_bbcsid, F_bbcsid] = bpdndf_dynamics_learning('BBC', [], [], inf_opts);              % Run the dictionary learning algorithm

%% BBC standard bpdn-df w/ 1 DL

inf_opts.special         = 'bpdn-df';
[D_1, F_1] = bpdndf_dynamics_learning('BBC', [], [], inf_opts);              % Run the dictionary learning algorithm

%% BBC standard bpdn-df w/ 6 DL

inf_opts.nF = 6;
inf_opts.special         = 'bpdn-df';
[D_6, F_6] = bpdndf_dynamics_learning('BBC', [], [], inf_opts);              % Run the dictionary learning algorithm

%% BBC standard bpdn-df w/ 12 DL

inf_opts.nF = 12;
inf_opts.special         = 'bpdn-df';
[D_12, F_12] = bpdndf_dynamics_learning('BBC', [], [], inf_opts);              % Run the dictionary learning algorithm

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
