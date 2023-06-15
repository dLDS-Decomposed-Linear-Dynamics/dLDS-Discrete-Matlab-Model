%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dictionary_learn_script.m 
% This scrips it an example script for basic sparse dictionary learning
% based on the Olshausen and Field, 1997 paper 'Sparse Coding with
% an Overcomplete Basis Set: a Strategy Employed by V1?'
% 
% Deviations from the model presented are: 
%   - use of L1 regularized BPDN solvers as well as options to use
%     hard-sparse MP-type solvers for the inference step. This allows
%     faster and more accurate solvutions for the sparse coefficients.
%   - Normalization of the dictionary elements to have unit norm at each
%     step rather than normalizing the variance of the coefficients. This
%     is the default. Also included is a method that normalizes the
%     dictionary elements by a Forbeneous norm.
%   
% Other optimizations:
% 
%   - Parallel for loops are used in the inference step to speed up
%     run-time. If the script is to be run on a cluster, ust the
%     createMatlabPoolJob() function to allow for this speedup. 
% 
% 
% Last Updated 6/3/2010 - Adam Charles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Options

% Set multithreading to Ncores - 1 threads
% This speeds up matrix multiplication etc.
maxNumCompThreads(8);
if matlabpool('size') ~= 0
    matlabpool CLOSE
end
matlabpool('open', 1)

% Random number seed
RandStream.setDefaultStream (RandStream('mt19937ar','seed',sum(100*clock)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options

% Save Options
opts.save_name = 'dict_learn_test.mat';  % Basis Save Name
opts.save_every = 500;              % How many iterations to save after

% Algorithm options
opts.data_type = 'square';          % Type of data (vecror, square or cube)
opts.sparse_type = 'l1ls';     % Choose to use l1ls for sparsification OMP_qr
                                    % Options: l1ls, l1ls_nneg, sparsenet, MP
opts.grad_type = 'norm';            % Choose weather to include the Forb norm in E(a,D)
opts.nneg_dict = 0;                 % Make the basis values all positive

% Dictionary options
opts.n_elem = 64;                   % Number of dictionary elements
opts.bl_size = 8;                   % Block Size (nxn)
opts.dep_size = 1;                  % Depth size for 'cube' data

% Iteration numbers
opts.iters = 1000;                  % Number of learning iterations
opts.in_iter = 75;                 % Number of internal iterations
opts.GD_iters = 1;                  % Basis Gradient Descent iterations

% Specify Parameters
opts.step_size = 5;                 % Initial Step Size for Gradient Descent
opts.decay = 0.9995;                % Step size decay factor
opts.lambda = 0.6;                  % Lambda Value for Sparsity
opts.lambda2 = 0.02;                % Forbenious norm lambda
opts.tol = 0.001;                   % Sparsification Tolerance
opts.h_sparse = 6;                  % Number of basis elements for hard-sparse

% Plotting Options
opts.bshow = 1;                     % Number of iterations to show basis
opts.disp_size = [8, 8];            % Basis display dimensions

% Data normalization options
opts.ssim_flag = 0;                 % Reduce variance for high variance blocks
opts.std_min = 0.10;                % Minimum standard deviation

%% Load Training Data

% Select data type. ''image'' or ''hyperspectral'' included here. Change to
% use other data types by simply loading the data into a matric called
% 'data_obj'. Check help file in learn_dictionary.m for data_type specific
% formats.

data_class = 'image'; 

if ~exist('data_obj', 'var')
    if strcmp(data_class, 'image')
        % Load Bruno's prewhitened Image set
        fprintf('Loading Images...\n')
        load('IMAGES.mat')
        data_obj = IMAGES;
        clear IMAGES;
        fprintf('Images Loaded\n')
    elseif strcmp(data_class, 'hyperspectral')
        % Note that the hyperspectral preprocessing depends on the type of
        % data available. The function 'hyper_bip_read' loads .BIP
        % formatted data and is not included in this code package. The
        % function 'hyper_data_fix' is also not included because it was
        % written for a particular dataset. For information on these
        % functions or to request the functions, email Adam at
        % acharles6@gatech.edu.
        header_file = 'P20011018_southendBIP.hdr';
        filename = 'P20011018_southendBIP';
        list_file = 'P20011018_southendBIP.wvl';
        filename_full = 'Probe_20011018_SmithI_tafkaa_effort_remap.dat';
        [data_obj, hdr_struct, hyper_index] = hyper_bip_read(filename, header_file, list_file);
        % [data_obj, hdr_struct, hyper_index] = hyper_bip_read(filename_full, list_file);
        h_opts.spec_num = size(data_obj, 1);
        h_opts.plotting = 0;            % Plotting Option for HS functions
        
        % Hyperspectral Preprocessing
        data_obj = hyper_data_fix(data_obj, h_opts);
        fprintf('Fixed Hyperspectral Data...\n')
        fprintf('Maximum value is %d, minimum value is %d, and mean value is %d.\n',...
            max(max(data_obj)), min(min(data_obj)),...
            sum(sum(data_obj))/(size(data_obj, 1)*size(data_obj, 2)))

    else
        error('Unknown Data Type!! Choose ''image'' or ''hyperspectral''...')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Dictionary Elements
% To load an older dictionary, 

if strcmp(opts.data_type, 'vector')
    v_size = size(data_obj, 2);
elseif strcmp(opts.data_type, 'square')
    v_size = opts.bl_size.^2;
elseif strcmp(opts.data_type, 'cube')
    v_size = opts.dep_size*(opts.bl_size.^2);
else
    error('Unknown Data Type!! Choose ''vector'', ''square'' or ''cube''...')
end

dictionary_initial = initialize_dictionary(opts.n_elem, v_size, opts.nneg_dict);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Learning Algorithm

% Actually run the learning algorithm
% [dictionary_out] = learn_dictionary(data_obj, dictionary_initial, opts);
[dictionary_out] = learn_dictionary_spmd(data_obj, dictionary_initial, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Closing time

matlabpool CLOSE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
