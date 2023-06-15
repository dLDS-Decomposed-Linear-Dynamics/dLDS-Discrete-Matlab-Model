function varargout = rand_seq_create(varargin)

% varargout = rand_seq_create(varargin)
% 
% Creates a random set of dynamics in cell form
% 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin < 2
    error('need to input both signal parameters and noise parameters')
else
    sig_params = varargin{1};
    noise_params = varargin{2};
end
if (~isfield(sig_params,'N'))||(~isfield(sig_params,'S'))||(~isfield(sig_params,'T_s'))
    error('need to specify signal size, sparsity and number of time steps')
end

if ~isfield(sig_params,'nF')
    sig_params.nF = 1;
elseif isempty(sig_params.nF)
    sig_params.nF = 1;
else                                                                       % do nothing
end

if ~isfield(sig_params,'sF')
    sig_params.sF = 1;
elseif isempty(sig_params.sF)
    sig_params.sF = 1;
else                                                                       % do nothing
end

if nargin > 2
    F_cell = varargin{3};
else
    F_cell = rand_dyn_create(sig_params.N, sig_params.nF, 'perm');
end

if ~isfield(sig_params,'sig_var')
    sig_params.sig_var = 1;
elseif isempty(sig_params.sig_var)
    sig_params.sig_var = 1;
else                                                                       % do nothing
end

if ~isfield(sig_params,'offset')
    sig_params.offset = 1;
elseif isempty(sig_params.offset)
    sig_params.offset = 1;
else                                                                       % do nothing
end

if ~isfield(sig_params,'unit_b')
    sig_params.unit_b = false;
elseif isempty(sig_params.unit_b)
    sig_params.unit_b = true;
else                                                                       % do nothing
end

if ~isfield(noise_params,'dg_var')
    noise_params.dg_var = 1e-4;
elseif isempty(noise_params.dg_var)
    noise_params.dg_var = 1e-4;
else                                                                       % do nothing
end

if ~isfield(noise_params,'p')
    noise_params.p = 0.2*sig_params.S;
elseif isempty(noise_params.p)
    noise_params.p = 0.2*sig_params.S;
else                                                                       % do nothing
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data

x     = zeros(sig_params.N,sig_params.T_s);
b     = ones(numel(F_cell),sig_params.T_s);
x_tmp = zeros(sig_params.N,1);
x_tmp(randsample(sig_params.N,sig_params.S)) = sqrt(sig_params.sig_var)*...
                              randn(sig_params.S, 1)+ sig_params.offset*...
                                               rand_bern(sig_params.S,0.5);
x(:,1)  = x_tmp;
tmp_idx = 1:sig_params.N;

for kk = 2:sig_params.T_s
    [F,b(:, kk)] = mk_bilinear_F(F_cell,sig_params.sF,sig_params.unit_b);  % Generate the weights and dynamics at the current iteration
    x_tmp        = F*x(:, kk-1);                                           % Pure prediction
    x_tmp(x_tmp~=0) = x_tmp(x_tmp~=0) + ...
                         sqrt(noise_params.dg_var)*randn(sum(x_tmp~=0),1); % Add Gaussian noise at each iteration


    if noise_params.p ~= 0                                                 % Check if there should be a sparse innovations
        s_idx   = tmp_idx(x_tmp~=0);                                       % Get the support of the signal
        s_nidx  = tmp_idx(x_tmp==0);                                       % Get the complimant to the support of the signal
        n_del   = max(min(rand_posn(1,1,noise_params.p),sig_params.S),0);  % Pick a small number of support elements to change
        idx_del = s_idx(randsample(numel(s_idx), n_del));                  % Pick locations to delete (from the support)
        idx_mk  = s_nidx(randsample(numel(s_nidx), n_del));                % Pick locations to delete (from the support compliment)
        
        x_tmp(idx_del) = 0;                                                % Initialize temporary version of the signal to 0
        x_tmp(idx_mk)  = sqrt(sig_params.sig_var)*randn(numel(idx_mk),1)...
                          + sig_params.offset*rand_bern(numel(idx_mk),0.5);% Change the support as per the locations sampled above
    else
    end
    x(:, kk) = x_tmp;                                                      % Save the current signal
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse outputs

if nargout > 0
varargout{1} = x;                                                          % Output the signal sequence
end
if nargout > 1
varargout{2} = F_cell;                                                     % Output the true dyanmics dictionary
end
if nargout > 2
varargout{3} = b;                                                          % Output the sequence of dynamics used
end
if nargout > 3
    for kk = 4:nargout
        varargout{kk} = [];
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
