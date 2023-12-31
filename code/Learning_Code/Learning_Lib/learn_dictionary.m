function [dictionary_end, coef_vals] = learn_dictionary(varargin)

% function [dictionary_end] = learn_dictionary(varargin)
% 
% This function learns a dictionary optimal for sparse coding according to
% the unsupervised, statistical method from Olshausen & Field 1996: 
% 
% for num_iterations
%   x = argmin ||y-Ax||^2 + lambda*||x||_1
%           OR
%   x = argmin ||y-Ax||^2 + lambda*||x||_0
% 
%   Take a step in the negative gradient of A & normalize the new vectors
% end
% 
% Inputs:
%   data_obj        - Data in some format: matrix, 3D array, or cell array
%                     For the 'vector' data type, data_obj should be a
%                     matrix with each column as a sample daat point. For
%                     the 'square' options, data_obj is a 3D array where
%                     each slice in the 3rd dimension is an example 2D
%                     array od data (e.g. concatenating images). For the
%                     'cube' data_type, data_obj should be a cell array,
%                     where each element in the array is an example (color
%                     image or video sequence) of the same size. Different
%                     sizes can cause errors when extracting the data.
%   initial_dict    - The starting dictionary for the algorithm. If left
%                     empty (using '[]'), learn_dictionary will initialize
%                     a dictionary with the appropriate size for the data
%                     type specified.
%   infer_handle    - The handle for a wrapper that calls an inference
%                     function. For more information see gen_multi_infer.m
%   opts            - Struct with various options, including:
%       .data_type  - Type of data ('vector', 'square' or 'cube')
%       .grad_type  - Type of cost function for gradient descent: 'norm' or
%                     'frob'.
%       .sparse_type- Type of inference to use. For full set of options,
%                     see multi_infer.m
%       .n_elem     - Number of dictionary elements
%       .in_iter    - Number of samples per iteration
%       .iters      - Number of iterations to run the algorithm for
%       .GD_iters   - Number of gradient descent steps per iteration
%       .verb       - 1 if verbose outputs are desired
%       .bl_size    - Block size for 'square' or 'cube' data types
%       .nneg_dict  - Choose nonegative values only for the dictionary
%       .step_size  - Initial step size for gradient descent
%       .decay      - Rate of decay for gradient descent step size
%       .lambda     - Lambda value for l1-regulated inference schemes
%       .lambda2    - Second lambda value for 'frob' energy function
%       .tol        - Tolerance for inference schemes
%       .ssim_flag  - Choose whether to normalize the data pre-inference
%       .std_min    - Minimum standard deviation (use with .ssim_flag)
%       .save_every - Number of iterations to save after
%       .save_name  - Filename to save the dictionary in
%       .bshow      - 0 to not plot intermediary dictionaries. Else plot
%                     every .bshow iterations.
%       .disp_size  - Dimensions of figure to display the dictionary
%
% Outputs
%   dictionary_end  - Final dictionary after the set number of iterations 
%   
% Uses: 
%   [dictionary_end] = learn_dictionary(data_obj, initial_dict, infer_handle, opts)
%           Learns a dictionaty on data_obj starting with the initial
%           dictionary given and using the options specified by the struct
%           opts. The inference function is specified by it's wrapper given
%           by infer_handle.
%   [dictionary_end] = learn_dictionary(data_obj, initial_dict, opts)
%           Learns a dictionaty on data_obj starting with the initial
%           dictionary given and using the options specified by the struct
%           opts. The inference type is defaulted to the conjugate gradient
%           method (@cg_l2l1_wrapper).
%   [dictionary_end] = learn_dictionary(data_obj, [], opts)
%   [dictionary_end] = learn_dictionary(data_obj, opts)
%           Learns a dictionary from the data in data_obj using the options
%           in opts. The initial dictionary is initialized randomly using
%           the function initialize_dictionary.m. The inference type is
%           defaulted to the conjugate gradient method (@cg_l2l1_wrapper).
% 
% DATA SETS TESTED ON:
% 
% There are two types of basis functions that have been tested with this
% algorithm. These are image patches and hyperspectral images.
% 
%   IMAGE PATCHES ('square')
% The data is input as an MxNxL data cube consisting of L images, each of
% size MxN. The images are assumed to be pre-whitened. 
% 
%   HYPERSPECTRAL DATA ('vector')
% The data is put in as a MxN (2D) matrix. For the basis learning, the
% image position in the image is irrelivent (no hierarchical model). The
% data should be 0 < data_{i,j} < 1. 
% 
% Last Updated 6/10/2010 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing - Set Defaults on the options

% Check the number of inputs and parse accordingly
if nargin == 2
    data_obj = varargin{1};
    opts = varargin{2};
    initial_dict = [];
    infer_handle = @l1ls_wrapper;
elseif nargin == 3
    data_obj = varargin{1};
    initial_dict = varargin{2};
    opts = varargin{3};
    infer_handle = @l1ls_wrapper;
elseif nargin == 4
    data_obj = varargin{1};
    initial_dict = varargin{2};
    infer_handle = varargin{3};
    opts = varargin{4};
end

% OPTIONS: Make sure that the correct options are set and that all
% necessary variables are available or defaulted.
if ~isfield(opts, 'verb')
    opts.verb = 1;                                                         % Default to no verbose output
end

if (~isfield(opts, 'data_type')) || (~isfield(opts, 'n_elem'))
    error('Missing vital options! You need at a minimum to include: data_type, n_elem.')
end

if (~isfield(opts, 'bl_size')) && ...
        (strcmp(opts.data_type, 'square') || strcmp(opts.data_type, 'cube'))
    error('Missing vital option! You need opts.bl_size for square or cube data!')
elseif (~isfield(opts, 'bl_size')) && strcmp(opts.data_type, 'vector')
    opts.bl_size = 1;
else
    if numel(opts.bl_size) == 1
        if strcmp(opts.data_type, 'square')
            opts.bl_size = opts.bl_size*ones(1,2);
        elseif strcmp(opts.data_type, 'cube')
            opts.bl_size = opts.bl_size*ones(1,3);
        else                                                               % Don't go here
        end
    elseif numel(opts.bl_size) == 2
        if strcmp(opts.data_type, 'square')                                % Everything is OK
        elseif strcmp(opts.data_type, 'cube')
            error('you gave 2 dimensions in bl_size for a 3D thing!')
        else                                                               % Don't go here
        end
    elseif numel(opts.bl_size) == 3
        if strcmp(opts.data_type, 'square')
            error('you gave 3 dimensions in bl_size for a 2D thing!')
        elseif strcmp(opts.data_type, 'cube')                              % Everything is OK
        else                                                               % Don't go here
        end
    else
        error('Too many dimensions in bl_size!')
    end
end
if ~isfield(opts, 'grad_type')
    opts.grad_type = 'norm';
end

if ~isfield(opts, 'in_iter')
    opts.in_iter = 200;                                                    % Default 200 samples per iteration
end

if ~isfield(opts, 'iters')
    opts.iters = 2000;                                                     % Default 2000 total iterations
end

if ~isfield(opts, 'ssim_flag')
    opts.ssim_flag = 0;                                                    % Default to no normalization between samples
end

if (~isfield(opts, 'std_min'))
    opts.std_min = 0.1;                                                    % Default to min. sample standard deviation of 0.1
    warning('Inputs:UnspecifiedParam',  ...
        ['Min sample STD not set by user!! Using STDmin = 0.1. ', ... 
        'This is probably bad! ', ...
        'Ctrl-C and restarting with a specified value is recommended.'])
end

if ~isfield(opts, 'tol')
    opts.tol = 0.001;                                                      % Default to tol = 0.001 (default tol for most optimizers)
end

if ~isfield(opts, 'GD_iters')
    opts.GD_iters = 1;                                                     % Default to one GD step per iteration
end

if ~isfield(opts, 'bshow')
    opts.bshow = 0;                                                        % Default to no plotting
end

if (~isfield(opts, 'disp_size')) && (opts.bshow >= 1)
    error('Missing display option, Must know dimensions of dictionary plotting!')
end

if ~isfield(opts, 'nneg_dict')
    opts.nneg_dict = 0;                                                    % Default to not having negativity constraints
end

if ~isfield(opts, 'lambda')
    opts.lambda = 0.1;                                                     % Default to lambda = 0.1
    warning('Inputs:UnspecifiedParam',  ...
        ['Lambda not set by user!! Using lambda = 0.1. ', ... 
        'This is probably bad for a non-image application! ', ...
        'Ctrl-C and restarting with a specified value is recommended.'])
end

if ~isfield(opts, 'step_size')
    opts.step_size = 3;                                                    % Default initial step size to 3
    warning('Inputs:UnspecifiedParam',  ...
        ['Step Size not set by user!! Using step size = 3. ', ... 
        'This is probably bad for a non-image application! ', ...
        'Ctrl-C and restarting with a specified value is recommended.'])
end

if ~isfield(opts, 'decay')
    opts.decay = 0.9995;                                                   % Default step size decay to 0.9995
    warning('Inputs:UnspecifiedParam',  ...
        ['Step Size decay not set by user!! Using decay of 0.9995. ', ... 
        'This is probably bad for a non-image application! ', ...
        'Ctrl-C and restarting with a specified value is recommended.'])
end

if (~isfield(opts, 'lambda2')) && (strcmp(opts.grad_type, 'frob'))
    opts.lambda2 = 0.02;                                                   % Default to lambda2 = 0.02
    warning('Inputs:UnspecifiedParam',  ...
        ['Lambda not set by user!! Using lambda2 = 0.02. ', ... 
        'This is probably bad for a non-image application! ', ...
        'Ctrl-C and restarting with a specified value is recommended.'])
end

if ~isfield(opts, 'save_every')
    opts.save_every = floor(opts.iters/10);                                % Default save every 10% iters
end

if ~isfield(opts, 'div_samp')
    opts.div_samp = 0;                                                     % Default to no resampling of similar exemplars in the SGD
end

if ~isfield(opts, 'div_thresh')
    opts.div_thresh = 0.85;                                                % Default threshold for sampling exemplars is 0.85
end

if ~isfield(opts, 'save_name')
    date_str = date;
    opts.save_name = [date_str(8:end), date_str(3:7), date_str(1:2), ...
        'Dictionary_' num2str(opts.n_elem), 'Elems_', num2str(opts.lambda), ...
        'lambda.mat'];
    fprintf('Save name not specified, saving as %s...\n', opts.save_name)
end

if ~isfield(opts, 'sparse_type')
    opts.sparse_type = 'l1ls';                                             % Default to l1ls for inference
end

% INITIAL DICTIONARY: Make sure that the dictionary is initialized. Only
% initializes randomly. For specific initial dictionary, load it prior to
% running and include as an input.
if isempty(initial_dict)                                                   % If no initialization given, set the initial dictionary here.
    if strcmp(opts.data_type, 'vector')
        v_size = size(data_obj, 1);
    elseif strcmp(opts.data_type, 'square')
        v_size = prod(opts.bl_size);
    elseif strcmp(opts.data_type, 'cube')
        v_size = prod(opts.bl_size);
    else
        error('Unknown Data Type!! Choose ''vector'', ''square'' or ''cube''...')
    end
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock))); % Random number seed
    initial_dict = initialize_dictionary(opts.n_elem, v_size, opts.nneg_dict); % Set the initial dictionary
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error Checking

if strcmp(opts.data_type, 'vector')
    if size(data_obj, 1) ~= size(initial_dict, 1)
        error('Dimension mismatch between data and initial dictionary size!')
    end
    if ~isnumeric(data_obj)
        error('data_obj must be in matrix form for ''vector'' type data!')
    end
    if numel(size(data_obj)) ~= 2
        error('data_obj must be in matrix form for ''vector'' type data!')
    end
elseif strcmp(opts.data_type, 'square')
    if prod(opts.bl_size) ~= size(initial_dict, 1)
        error('Dimension mismatch between opts.bl_size and initial dictionary size!')
    end
    if ~isnumeric(data_obj)
        error('data_obj must be in matrix form for ''square'' type data!')
    end
    if numel(size(data_obj)) ~= 3
        error('data_obj must be in matrix form for ''square'' type data!')
    end
    if (size(data_obj, 1) < opts.bl_size(1)) || (size(data_obj, 2) < opts.bl_size(2))
        error('The dimensions of the data objects is too small!')
    end
elseif strcmp(opts.data_type, 'cube')
    if opts.prod(opts.bl_size) ~= size(initial_dict, 1)
        error('Dimension mismatch between opts.bl_size and initial dictionary size!')
    end
    if ~iscell(data_obj)
        error('data_obj must be in cell form for ''cube'' type data!')
    end
    for ii = 1:numel(data_obj)
        if (size(data_obj, 1) < opts.bl_size(1)) || (size(data_obj, 2) < opts.bl_size(2)) ...
                || (size(data_obj, 3) < opts.bl_size(3))
            error('The dimensions for one or more of the data objects is too small!')
        end
    end
else
    error('Invalid data type! choose ''vector'', ''square'', or ''cube''!')
end

if opts.n_elem ~= size(initial_dict, 2)
    error('Dimension mismatch between opts.n_elem and initial dictionary size!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations and Dimention Extraction


dictionary_n = initial_dict;                                               % Initialize Basis
iter_num     = 0;                                                          % Iteration counter initialization
step_s       = opts.step_size;                                             % Initialize step size
% x_im = zeros(size(dictionary_n,1), opts.in_iter);                          % Initialize matricies that will be populated during the actual learning
% std_min_vec = opts.std_min;                                                % Create variable so that MATLAB's parallel toolbox doesn't crash
RandStream.setGlobalStream (RandStream('mt19937ar','seed',sum(100*clock)));% Random number seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Algorithm

fprintf('Educating your basis...\n')
basic_cell.options = opts;

tic
while iter_num <= opts.iters
    try
        iter_num = iter_num + 1;                                           % Update iteration number
        
        %% Get Training Data. 
        % Some parts can be run in parallel with appropriate slicing.
        
        x_im = sample_exemplars(data_obj, opts);
        
        %% Interior loop: find sparse coefficients 
        % There are a number of packages to accomplish this. This is one
        % block of code that can be MUCH improved with parallalization.

        coef_vals = gen_multi_infer(dictionary_n, x_im, infer_handle, opts);

        %% Minimize the energy w.r.t. the dictionary using gradient descent
        
       	dictionary_n = dictionary_update(x_im, dictionary_n, coef_vals, step_s, opts);

        if opts.verb == 1
            % Spit out info
            im_snr = mean(sum((x_im - dictionary_n*coef_vals).^2, 1)./...
                                                          sum(x_im.^2, 1));
            sparse_num = mean(sum(abs(coef_vals) <= 1e-7));
            fprintf('Iteration %d, SNR is %d, sparsity is %d, step size is %d, max coeffs is %d \n',...
                iter_num, im_snr, sparse_num, step_s, max(max(coef_vals)))
        end
        
        step_s = step_s*opts.decay;                                        % Update the step size
        
        %% Plot Outputs: Look at basis elements every opts.bshow iterations
        if opts.bshow ~= 0 && (mod(iter_num, opts.bshow) == 0)
            if strcmp(opts.data_type, 'vector')
                fprintf('Plotting...\n')
                % Plot the columns of dictionary_n
                x_vect = 1:size(dictionary_n, 1);
                dict_plot1d(dictionary_n, x_vect, opts.disp_size, 100)
                drawnow
            elseif strcmp(opts.data_type, 'square')
                % Make basis into image patches for display purposes
                I = basis2img2(dictionary_n, [opts.bl_size,opts.bl_size], opts.disp_size);
                % Display Image with title Etc.
                imagesc(I)
                title(sprintf('Learned Dictionary: Iteration %d, Image Patches',...
                    iter_num), 'FontSize', 16);
                drawnow
            elseif strcmp(opts.data_type, 'cube')
                if opts.bl_size(3) == 3
                    fprintf('Plotting as Color Image...\n')
                    I1 = basis2img2(dictionary_n(1:opts.bl_size^2, :), [opts.bl_size,opts.bl_size], opts.disp_size);
                    I2 = basis2img2(dictionary_n(opts.bl_size^2+1:2*opts.bl_size.^2, :),...
                        [opts.bl_size,opts.bl_size], opts.disp_size);
                    I3 = basis2img2(dictionary_n(2*opts.bl_size.^2+1:3*opts.bl_size.^2, :),...
                        [opts.bl_size,opts.bl_size], opts.disp_size);
                    I(:,:,1) = I1(:,:,1);
                    I(:,:,2) = I2(:,:,1);
                    I(:,:,3) = I3(:,:,1);
                    figure(1)
                    imshow(I, 'InitialMagnification', 300)
                    drawnow
                else
                    warning('Plotting:DataType', 'This is not a color image, so plotting is futile!')
                end
            else
                error('Unknown Data Type!!')
            end
        else
            % Do nothing
        end

        
        %% Save data... 
        if mod(iter_num, opts.save_every) == 0
            fprintf('Saving progress...\n')
            basic_cell.dictionary = dictionary_n;
            basic_cell.iter = iter_num;
            eval(sprintf('save %s basic_cell;', opts.save_name));
        else
            % Do nothing
        end
    catch ME
        fprintf('Saving last dictionary before error...\n')
        basic_cell.dictionary = dictionary_n;
        basic_cell.iter = iter_num;
        eval(sprintf('save %s basic_cell;', opts.save_name));
        
        fprintf(ME.message)
        fprintf('The program failed. Your dictionary at the last iteration was saved.')
        rethrow(ME)
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output Stuff and wrapping up
tot_time = toc;
        
% Save everything
fprintf('Saving final state...\n')
basic_cell.dictionary = dictionary_n;
basic_cell.iter = iter_num;
eval(sprintf('save %s basic_cell;', opts.save_name));

fprintf('Total iteration time was %d\n', tot_time);

dictionary_end = dictionary_n;

if nargout == 2
    opts.in_iter = size(data_obj,2);
    coef_vals = gen_multi_infer(dictionary_end, data_obj, infer_handle, opts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
