
clear

sig_opts.N = 20;
sig_opts.S = 5;
sig_opts.sig_var = 0.05;
sig_opts.offset = 1;
sig_opts.T_s = 6;

noise_opts.p = 1;
noise_opts.dg_var = 1e-4;

[x, F] = rand_seq_create(sig_opts, noise_opts);


imagesc(x)