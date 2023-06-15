function varargout = Load_Data(varargin)


if nargout < 3
    error('not enough output arguments')
end

if nargin == 0
    error('not enough input arguments')
end

data_type = varargin{1};

if strcmp(data_type, 'yuv_vid') == 1
    % Load YUV video and simulate CS measurements
    [y, MEAS_FUN, x0] = VID_loadNsamp(varargin{2}, varargin{3});
elseif strcmp(data_type, 'mri') == 1
    % Load MRI and simulate CS measurements
    [y, MEAS_FUN, x0] = MRI_loadNsamp(varargin{2}, varargin{3});
elseif strcmp(data_type, 'track') == 1
    param_obj = varargin{2};
    if (~isfield(param_obj,'M'))||(~isfield(param_obj,'N'))||(~isfield(param_obj,'T_s'))||(~isfield(param_obj,'S'))
        error('Need to specify M, N, and T_s')
    end
    if ~isfield(param_obj,'n_paths')
        param_obj.n_paths = [];
    end
    if ~isfield(param_obj,'noise_var')
        param_obj.noise_var = [];
    end
    if ~isfield(param_obj,'one_step')
        param_obj.one_step = [];
    end

    for kk = 1:num_trials
        [x0, F] = gen_path_seq(param_obj.S, param_obj.T_s, sqrt(param_obj.N)*[1, 1],  param_obj.p, 1);
        [y, G] = take_gaussian_meas(x0, param_obj.M, param_obj.noise_var);
    end
    MEAS_FUN = cell(1, T_s);
    DYN_FUN = cell(1, T_s+1);

    for kk = 1:param_obj.T_s
        meas_func.Phi  = @(z) G(:,:,kk)*z(:);
        meas_func.Phit = @(z) (G(:,:,kk)')*z(:);
        MEAS_FUN{kk} = meas_func;
        
	f_dyn = @(z) F(:,:,kk)*z;
	DYN_FUN{kk+1} = f_dyn;
    end
elseif strcmp(data_type, 'com_chan') == 1
    % Generate communications channel data
    param_obj = varargin{2};
    if (~isfield(param_obj,'M'))||(~isfield(param_obj,'N'))||(~isfield(param_obj,'T_s'))
        error('Need to specify M, N, and T_s')
    end
    if ~isfield(param_obj,'n_paths')
        param_obj.n_paths = [];
    end
    if ~isfield(param_obj,'noise_var')
        param_obj.noise_var = [];
    end
    if ~isfield(param_obj,'one_step')
        param_obj.one_step = [];
    end

    [y, MEAS_FUN, x0] = make_com_chan_data(param_obj.M, param_obj.N, param_obj.T_s, ...
        param_obj.n_paths, param_obj.noise_var, param_obj.one_step);

    % DWTfunc.apply = @(q) q;
    % DWTfunc.invert = @(q) q;
else
    error('unknown data type to load')
end

varargout{1} = y;
varargout{2} = MEAS_FUN;
varargout{3} = x0;
if nargout > 3
    if ~exists('DYN_FUN')
        varargout{4} = [];
    else
        varargout{4} = DYN_FUN;
    end
end

