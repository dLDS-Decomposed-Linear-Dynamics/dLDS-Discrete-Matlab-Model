%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Script that tests dictionary learning for BPDNDF dynamics %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set parameters
block_size = 4;

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
X_ex = X_ex(:,301:end);

%%

inf_opts.M               = sig_opts.M;
inf_opts.N               = sig_opts.N;
inf_opts.lambda_val      = 0.01;
inf_opts.lambda_history  = 0.1;
inf_opts.lambda_b        = 0.001;
% inf_opts.lambda_history  = 0.90;
% inf_opts.lambda_b        = 0.01;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.N_ex            = 100;
inf_opts.nF              = 2;
inf_opts.max_iter        = 10e3;
inf_opts.max_iter2       = 100;
inf_opts.step_d          = 70;
inf_opts.step_f          = 10;
inf_opts.plot_option     = 1;

inf_opts.special         = 'bpdn-df';
[D, F] = bpdndf_dynamics_learning(10*X_ex, [], [], inf_opts);              % Run the dictionary learning algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%