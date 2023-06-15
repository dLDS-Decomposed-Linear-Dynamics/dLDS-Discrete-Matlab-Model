%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options


sig_opts.N = 64;
sig_opts.M = 64;
sig_opts.S = 10;
sig_opts.sig_var = 0.05;
sig_opts.sig_var = 1;
sig_opts.T_s = 20;
sig_opts.nF = 10;
sig_opts.sF = 2;

noise_opts.p = 0;
noise_opts.dg_var = 1e-4;

inf_opts.lambda_val = 0.01;
inf_opts.gamma_val = 0.01;
inf_opts.alpha_a = 0.2;
inf_opts.beta_a = 2;
inf_opts.xi_a = 1;
inf_opts.alpha_b = 0.2;
inf_opts.beta_b = 2;
inf_opts.tol = 1e-3;

N_ex = 5;

infer_hand = @rwl1df_bilinear_handle;

F_true = rand_dyn_create(sig_opts.N, sig_opts.nF, 'perm');
D_true = eye(sig_opts.N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

step_f = 0*20;
step_d = 0*10;
step_decay = 0.9995;

F = initialize_dynamics(sig_opts.nF, sig_opts.N);
D = initialize_dictionary(sig_opts.N, sig_opts.M);

% For the time being just test the inference!
F_init = F_true;
D_init = D_true;
% F_init = F;
% D_init = D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Learning

for n_iters = 1:2

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make some data

    X_ex = cell(1,N_ex);

    for kk = 1:N_ex
        [X_ex{kk}, ~, ~] = rand_seq_create(sig_opts, noise_opts, F_true);
	%X_ex{kk} = reshape(D_true*X_ex{kk},size(D,1),1,size(X_ex{kk},2)); 
	X_ex{kk} = D_true*X_ex{kk}; 
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
    
    if 10*floor(n_iters/10) == n_iters
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
    end
    
    % Display some statistics
    f_err = correlate_learned_dynamics(F, F_true,D);
    
    fprintf('Iteration %d complete. D error is %f and F error is %f (avg b is %f)\n', n_iters, ...
        sum((vec(D*(D')) - D_true(:)).^2), mean(f_err), mean(mean(cell2mat(B_cell))) )
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
