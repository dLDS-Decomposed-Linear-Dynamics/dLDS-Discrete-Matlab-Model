%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Force Start

fprintf('About to open...\n')
% [poolsize] = PaceParallelToolbox_r2014a_2(1);

% addpath(genpath('/nv/hp16/acharles6/VersionedCode/2013_Flearning/'));
addpath(genpath(fullfile(hm,r_dir,'dynamics_learning','code')))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options

d_v = clock;
plot_every = 50;
save_every = 50;
load_training_flag = true;

%% Options

sig_opts.N = 64;
sig_opts.M = 64;
sig_opts.S = 10;
sig_opts.sig_var = 0.05;
sig_opts.sig_var = 1;
sig_opts.T_s = 20;
sig_opts.nF = 16;
sig_opts.sF = 3;

noise_opts.p = 0;
noise_opts.dg_var = 1e-4;

inf_opts.lambda_val = 0.01;
inf_opts.lambda_history = 0.10;
inf_opts.lambda_b = 1e-4;
inf_opts.lambda_historyb = 0;
inf_opts.tol = 1e-3;
inf_opts.special = '';

N_ex = 40;

infer_hand = @bpdndf_bilinear_handle;
data_type = 'synth';

if strcmp(data_type, 'synth')
    F_true = rand_dyn_create(sig_opts.N, sig_opts.nF, 'perm');
    D_true = eye(sig_opts.N);
elseif strcmp(data_type, 'video')
    % Video, nothing to initialize
else
    % ???
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

step_f = 20;
step_d = 10;
step_decay = 0.99995;

F = initialize_dynamics(sig_opts.nF, sig_opts.N);
D = initialize_dictionary(sig_opts.N, sig_opts.M);

F_init = F;
D_init = D;

if strcmp(data_type, 'synth')
    % save_name = sprintf('/nv/hp16/acharles6/data/%4.0f%2.0f%2.0fFORCE_test_synth_nf%2.0f_sf%2.0f.mat',d_v(1),d_v(2),d_v(3),sig_opts.nF,sig_opts.sF);
    save_name = 'synth_rerun_FORCE_obs.mat';
elseif strcmp(data_type, 'video')
    save_name = sprintf('/nv/hp16/acharles6/data/%4.0f%2.0f%2.0fFORCE_test_bbc_nf%2.0f_sf%2.0f.mat',d_v(1),d_v(2),d_v(3),sig_opts.nF,sig_opts.sF);
else
end

if load_training_flag
    load(save_name)
    start_iter = n_iters;
else
    start_iter = 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Learning

for n_iters = start_iter:20000

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make some data

    X_ex = cell(1,N_ex);

    if strcmp(data_type, 'synth')
        for kk = 1:N_ex
            [X_ex{kk}, ~, ~] = rand_seq_create(sig_opts, noise_opts, F_true);
	    %X_ex{kk} = reshape(D_true*X_ex{kk},size(D,1),1,size(X_ex{kk},2)); 
	    X_ex{kk} = D_true*X_ex{kk}; 
        end
    elseif strcmp(data_type, 'video')
        % Get random videos
	for kk = 1:N_ex
	    X_ex{kk} = rand_bbc_video([sqrt(sig_opts.N), sqrt(sig_opts.N), sig_opts.T_s]);
	end
    else
        error('Unknown data type!')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Inference

    [A_cell,B_cell] = parallel_bilinear_dynamic_inference(X_ex, D, F, infer_hand, inf_opts);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Update step

    D = dictionary_update(cell2mat(X_ex), D, cell2mat(A_cell), step_d);
    if sig_opts.T_s > 1 
        F = dynamics_update(F, A_cell, step_f,B_cell);
    end
    step_d = step_d*step_decay;
    step_f = step_f*step_decay;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plotting?
    
    if strcmp(data_type, 'synth')
    if (plot_every*floor(n_iters/plot_every) == n_iters)&&(plot_every > 0)
        %%
        figure(1)
        subplot(2,2,1), imagesc(basis2img2(D, sqrt(sig_opts.N)*[1,1], sqrt(sig_opts.N)*[1,1]))
        title('Learned Dictionary', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        subplot(2,2,2), imagesc((D)*F{1}*(D'))
        title('Learned Dynamics', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        figure(1)
        subplot(2,2,3), imagesc(basis2img2(D_true, sqrt(sig_opts.N)*[1,1], sqrt(sig_opts.N)*[1,1]))
        title('True Dictionary', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        subplot(2,2,4), imagesc(F_true{1})
        title('True Dynamics', 'FontSize', 20)
        colormap gray
        axis image
        axis off
    else
        % No Plot
    end
    if (save_every*floor(n_iters/save_every) == n_iters)&&(save_every > 0)
        save(save_name,'F','D','F_init','D_init','F_true','D_true','sig_opts','noise_opts','inf_opts','n_iters','N_ex','infer_hand','data_type','d_v');
    else
        % No Save
    end
    elseif strcmp(data_type, 'video')
    if (plot_every*floor(n_iters/plot_every) == n_iters)&&(plot_every > 0)
        %%
        figure(1)
        subplot(1,2,1), imagesc(basis2img2(D, sqrt(sig_opts.N)*[1,1], sqrt(sig_opts.N)*[1,1]))
        title('Learned Dictionary', 'FontSize', 20)
        colormap gray
        axis image
        axis off
	for mm = 1:numel(F)
	    Fmat(:, mm) = reshape(F{mm}, [], 1);
	end
        subplot(1,2,2), imagesc(basis2img2(Fmat, size(F{1}), [4,5]))
        title('Learned Dynamics', 'FontSize', 20)
        colormap gray
        axis image
        axis off
    else
        % No Plot
    end
    if (save_every*floor(n_iters/save_every) == n_iters)&&(save_every > 0)
        save(save_name,'F','D','F_init','D_init','sig_opts','noise_opts','inf_opts','n_iters','N_ex','infer_hand','data_type','d_v');
    else
        % No Save
    end
    else
        error('Unknown data type!')
    end

    if strcmp(data_type, 'synth')
        % Display some statistics
        f_err = correlate_learned_dynamics(F, F_true,D);
    
        fprintf('Iteration %d complete. D error is %f and F error is %f (avg b is %f)\n', n_iters, ...
            sum((vec(D*(D')) - D_true(:)).^2), mean(f_err), mean(mean(cell2mat(B_cell))) )
    elseif strcmp(data_type, 'video')
        fprintf('Iteration %d complete. (avg b is %f)\n', n_iters, mean(mean(cell2mat(B_cell))))
    else
        error('Unknown data type!')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final Save
if strcmp(data_type, 'synth')
    save(save_name,'F','D','F_init','D_init','F_true','D_true','sig_opts','noise_opts','inf_opts','n_iters','N_ex','infer_hand','data_type','d_v');
elseif strcmp(data_type, 'video')
    save(save_name,'F','D','F_init','D_init','sig_opts','noise_opts','inf_opts','n_iters','N_ex','infer_hand','data_type','d_v');
else
        error('Unknown data type!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Force

delete(gcp('nocreate'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
