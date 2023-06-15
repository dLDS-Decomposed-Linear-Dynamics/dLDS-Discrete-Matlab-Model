function [D, F] = bpdndf_dynamics_learning(data_obj, F_init, D_init, inf_opts)

% [D, F] = bpdndf_dynamics_learning(data_obj, F_init, D_init, inf_opts)
% 
% Function to learn both the representation and dynamics dictionaries for
% the BPDN-DF model.
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse Inputs

[data_type, dataShape, sig_opts] = checkInputTypes(data_obj);              % Check the input type and set up some parameters accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter checking

inf_opts = check_inf_params(inf_opts);                                     % Check all the basic parameters
inf_opts = checkSizes(inf_opts, D_init, F_init);                           % Check the sizes of all the initializations
[D_true, F_true] = generateOptionalGT(data_type, sig_opts); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

if isempty(F_init);   F = initialize_dynamics(inf_opts.nF, inf_opts.N);    % If necessary, initialize the dynamics
else;                 F = F_init;                                          % Otherwise just pass through
end
if isempty(D_init);   D = initialize_dictionary(inf_opts.N, inf_opts.M);   % If necessary, initialize the dictionary
else;                 D = D_init;                                          % Otherwise just pass through
end

clear D_init F_init                                                        % Do some house cleaning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Learning

for n_iters = 1:inf_opts.max_iters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make some data
    
    X_ex = sampleSomeDataSeqs(inf_opts, data_type, sig_opts, data_obj, ...
                                                           D_true, F_true);% Select some data for this learning iteration
                                                       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Inference

    [A_cell,B_cell] = parallel_bilinear_dynamic_inference(X_ex, D, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Update step
    
    D_old = D;  F_old = F;                                                 % Retain a copy of the representational & dynamics dictionary for comparison purposes
    if (inf_opts.D_update == true) % (~strcmp(inf_opts.special, 'noobs'))||   % ASC added 7/25
        D = dictionary_update(cell2mat(X_ex), D, cell2mat(A_cell), ...
                                                         inf_opts.step_d); % Update static dictionary
    end  % ASC added 7/25
    if (inf_opts.T_s > 1 && isempty(inf_opts.F_update)) || (inf_opts.F_update == true)
        F = dynamics_update(F, A_cell, inf_opts.step_f, B_cell, inf_opts.lambda_f);           % Update dynamics dictionary
    end
    inf_opts.step_d = inf_opts.step_d*inf_opts.step_decay;                 % Reduce the step size for the static dictionary
    inf_opts.step_f = inf_opts.step_f*inf_opts.step_decay;                 % Reduce the step size for the dynamics dictionary
    
    inf_opts.lambda_f = inf_opts.lambda_f*inf_opts.lambda_f_decay;         % Reduce the regularization for the dynamics dictionary
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Display some statistics
    outputLearningStats(data_type,D,F,inf_opts,n_iters,D_true,F_true,D_old,F_old,B_cell,dataShape);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRA FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [data_type, dataShape, sig_opts] = checkInputTypes(data_obj)

sig_opts = [];                                                             % Initialize option struct
if isstruct(data_obj)
    sig_opts  = data_obj;
    data_type = 'synth';
    dataShape = 'image';
elseif ischar(data_obj)
    data_type = data_obj;
    if strcmp(data_obj,'synth')
        sig_opts = default_synth_params();                                 % Get default parameters
    end
    dataShape = 'image';
elseif iscell(data_obj)
    data_type = 'datacell';
    if ismatrix(data_obj{1})
        dataShape = 'vector';
    elseif ndims(data_obj{1}) == 3
        dataShape = 'image';
    else
        error('Incompatiable number of dimensions in the data array!')
    end
elseif isnumeric(data_obj)
    data_type = 'datamatrix';
    if ismatrix(data_obj)
        dataShape = 'vector';
    elseif ndims(data_obj) == 3
        dataShape = 'image';
    else
        error('Incompatiable number of dimensions in the data array!')
    end
    
elseif isempty(data_obj)
    data_type = 'synth';
    dataShape = 'image';
    sig_opts  = default_synth_params();
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function inf_opts = checkSizes(inf_opts, D_init, F_init)

if (~isfield(inf_opts,'M'))||isempty(inf_opts.M)
    if ~isempty(D_init)
        inf_opts.M = size(D_init,1);
    elseif strcmp(data_type,'datamatrix')
        inf_opts.M = size(data_obj,1);
    elseif isfield(sig_opts,'M')
        inf_opts.M = sig_opts.M;
    else
        inf_opts.M = 144;
    end
end
if (~isfield(inf_opts,'N'))||isempty(inf_opts.N)
    if ~isempty(D_init)
        inf_opts.N = size(D_init,2);
    elseif isfield(sig_opts,'N')
        inf_opts.N = sig_opts.N;
    else
        inf_opts.N = 4*inf_opts.M;
    end
end
if (~isfield(inf_opts,'nF'))||isempty(inf_opts.nF)
    if ~isempty(F_init)
        inf_opts.nF = numel(F_init);
    elseif isfield(sig_opts,'nF')
        inf_opts.nF = sig_opts.nF;
    else
        inf_opts.nF = 25;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function X_ex = sampleSomeDataSeqs(inf_opts, data_type, sig_opts, data_obj, D_true, F_true)

X_ex = cell(1,inf_opts.N_ex);
if strcmp(data_type, 'synth')
    for kk = 1:inf_opts.N_ex
        [X_ex{kk}, ~, ~] = rand_seq_create(sig_opts, ...
                                         sig_opts.noise_opts, F_true); % Synthesize a new set of data examples
        X_ex{kk} = D_true*X_ex{kk}; 
    end
elseif strcmp(data_type, 'BBC')
    for kk = 1:inf_opts.N_ex
        X_ex{kk} = rand_bbc_video([sqrt(inf_opts.M), ... 
                     sqrt(inf_opts.M), inf_opts.T_s], inf_opts.dsamp); % Sample some videos from the BBC dataset
    end
elseif strcmp(data_type, 'datacell')||strcmp(data_type, 'datamatrix')
    X_ex = sample_dynamic_exemplars(data_obj, inf_opts);
else
    error('Unknown data type')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [D_true, F_true] = generateOptionalGT(data_type, sig_opts)
if strcmp(data_type, 'synth')
    F_true = rand_dyn_create(sig_opts.M, sig_opts.nF, 'perm');             % Create a dictionary of dynamics
    D_true = eye(sig_opts.M);                                              % Make a simple test case of a dictionary of points
else     
    F_true = []; D_true = [];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function outputLearningStats(data_type,D,F,inf_opts,n_iters,D_true,F_true,D_old,F_old,B_cell,dataShape)

if strcmp(data_type, 'synth')
    tmp_plotting_fun(data_type,D,F,inf_opts,n_iters,D_true,F_true)
    f_err = correlate_learned_dynamics(F, F_true,D, dataShape);

    fprintf('Iteration %d complete. D error is %f and F error is %f (avg b is %f)\n', n_iters, ...
        sum((vec(D*(D')) - D_true(:)).^2), mean(f_err), mean(mean(cell2mat(B_cell))) )
elseif strcmp(data_type, 'BBC')||strcmp(data_type, 'datamatrix')||strcmp(data_type, 'datacell')
    tmp_plotting_fun(data_type,D,F,inf_opts,n_iters,dataShape)
    f_err = computeDynErr(F,F_old);
    fprintf('Iteration %d complete. Mean dD: %f, Mean dF: %f, Max b: %f, Median b use: %f\n', n_iters, ...
        mean(sum((D - D_old).^2)./sum(D_old.^2)), f_err, max(max(abs(cell2mat(B_cell)))), ...
        median(sum(abs(cell2mat(B_cell))>1e-3)) )
else
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function f_err = computeDynErr(F,F_old)
f_err = 0;
for mm = 1:numel(F)
    f_err = f_err + (1/numel(F))*sum(sum((F_old{mm} - F{mm}).^2))/sum(sum(F{mm}.^2));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%