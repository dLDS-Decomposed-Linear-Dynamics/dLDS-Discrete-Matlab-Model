%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options

block_size = 12;

sig_opts.N = 4*block_size^2;
sig_opts.M = block_size^2;
sig_opts.S = 10;
sig_opts.sig_var = 0.05;
sig_opts.sig_var = 1;
sig_opts.T_s = 20;
sig_opts.nF = 25;
sig_opts.sF = 2;
sig_opts.dsamp = 4;

noise_opts.p = 0;
noise_opts.dg_var = 1e-4;

inf_opts.lambda_val = 0.8;
inf_opts.lambda_history = 0.90;
inf_opts.lambda_b = 0.02;
inf_opts.lambda_historyb = 0;
inf_opts.tol = 1e-3;

plot_option = 100;
N_ex = 10;

infer_hand = @bpdndf_bilinear_handle;
data_type = 'video';

if strcmp(data_type, 'synth')
    F_true = rand_dyn_create(sig_opts.M, sig_opts.nF, 'perm');
    D_true = eye(sig_opts.M);
elseif strcmp(data_type, 'video')
    % Video, nothing to initialize
else
    % ???
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

step_f = 20;
step_d = 10;
step_decay = 0.9998;

F = initialize_dynamics(sig_opts.nF, sig_opts.N);
D = initialize_dictionary(sig_opts.N, sig_opts.M);

F_init = F;
D_init = D;

Fmat = zeros(numel(F{1}), numel(F));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Learning

for n_iters = 1:30000

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make some data
    
    X_ex = cell(1,N_ex);

    if strcmp(data_type, 'synth')
        for kk = 1:N_ex
            [X_ex{kk}, ~, ~] = rand_seq_create(sig_opts, noise_opts, F_true);
	        % X_ex{kk} = reshape(D_true*X_ex{kk},size(D,1),1,size(X_ex{kk},2)); 
	        X_ex{kk} = D_true*X_ex{kk}; 
        end
    elseif strcmp(data_type, 'video')
        % Get random videos
	for kk = 1:N_ex
	    X_ex{kk} = rand_bbc_video([sqrt(sig_opts.M), sqrt(sig_opts.M), sig_opts.T_s], sig_opts.dsamp);
	end
    else
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Inference

    [A_cell,B_cell] = parallel_bilinear_dynamic_inference(X_ex, D, F, infer_hand, inf_opts);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Update step
    D_old = D;
    F_old = F;
    D = dictionary_update(cell2mat(X_ex), D, cell2mat(A_cell), step_d);
    if sig_opts.T_s > 1 
        F = dynamics_update(F, A_cell, step_f,B_cell);
    end
    step_d = step_d*step_decay;
    step_f = step_f*step_decay;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plotting?
    
    if strcmp(data_type, 'synth')
    if plot_option*floor(n_iters/plot_option) == n_iters
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
        drawnow
    end
    elseif strcmp(data_type, 'video')
    if plot_option*floor(n_iters/plot_option) == n_iters
        %%
        figure(1)
        subplot(1,2,1), imagesc(basis2img2(D, sqrt(sig_opts.M)*[1,1], sqrt(sig_opts.N)*[1,1]))
        title('Learned Dictionary', 'FontSize', 20)
        colormap gray
        axis image
        axis off
	    for mm = 1:numel(F)
	            Fmat(:, mm) = reshape(F{mm}, [], 1);
	    end
        subplot(1,2,2), imagesc(basis2img2(Fmat, size(F{1}), [5,5]))
        title('Learned Dynamics', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        drawnow
    end
    else
    end
    %% Display some statistics
    if strcmp(data_type, 'synth')
        f_err = correlate_learned_dynamics(F, F_true,D);
        
        fprintf('Iteration %d complete. D error is %f and F error is %f (avg b is %f)\n', n_iters, ...
            sum((vec(D*(D')) - D_true(:)).^2), mean(f_err), mean(mean(cell2mat(B_cell))) )
    elseif strcmp(data_type, 'video')
        f_err = 0;
        for mm = 1:numel(F)
            f_err = f_err + (1/numel(F))*sum(sum((F_old{mm} - F{mm}).^2))/sum(sum(F{mm}.^2));
        end
        fprintf('Iteration %d complete. Mean dD: %f, Mean dF: %f, Max b: %f, Median b use: %f\n', n_iters, ...
            mean(sum((D - D_old).^2)./sum(D_old.^2)), f_err, max(max(abs(cell2mat(B_cell)))), ...
            median(sum(abs(cell2mat(B_cell))>1e-3)) )
    else
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
