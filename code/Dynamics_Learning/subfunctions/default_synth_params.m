function sig_opts = default_synth_params()

sig_opts.block_size        = 12;
sig_opts.N                 = 4*sig_opts.block_size^2;
sig_opts.M                 = sig_opts.block_size^2;
sig_opts.S                 = 10;
sig_opts.sig_var           = 0.05;
sig_opts.sig_var           = 1;
sig_opts.T_s               = 20;
sig_opts.nF                = 25;
sig_opts.sF                = 2;
sig_opts.dsamp             = 4;
sig_opts.noise_opts.p      = 0;
sig_opts.noise_opts.dg_var = 1e-4;

end