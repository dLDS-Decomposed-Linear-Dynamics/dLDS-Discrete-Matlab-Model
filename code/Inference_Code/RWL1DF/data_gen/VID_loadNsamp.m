function varargout = VID_loadNsamp(varargin)

filename = varargin{1};
param_obj = varargin{2};

if (~isfield(param_obj,'samp_factor'))||(~isfield(param_obj,'T_s'))
    error('Need to specify samp_factor and T_s')
end

samp_factor = param_obj.samp_factor;
T_s = param_obj.T_s;

if ~isfield(param_obj,'DCT_fact')
    DCT_fact = 0;
else
    DCT_fact = param_obj.DCT_fact;
    if isempty(DCT_fact)
        DCT_fact = 0;
    end
end
if ~isfield(param_obj,'noise_var')
    noise_var = 0.0001;
else
    noise_var = param_obj.noise_var;
    if isempty(noise_var)
        noise_var = 0.0001;
    end
end



%if ~isfield(param_obj,'num_trials')
%    num_trials = 1;
%else
%    num_trials = param_obj.num_trials;
%    if isempty(num_trials)
%        num_trials = 1;
%    end
%end
% normalize_opt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load video
% vid_t = yuv_import(filename, [176 144], T_s*num_trials);
vid_t = yuv_import(filename, [176 144], T_s);

% Extract the signals of interest
vid = cell(T_s);
for kk = 1:T_s
    % vid{kk} = vid_t{T_s*(trial_index-1)+kk};
    vid{kk} = vid_t{kk};
    vid{kk} = vid{kk}(1:128, 1:128);
end
clear vid_t

N = numel(vid{1});
M = ceil(samp_factor*N);

K_DCT = floor(0.5*(-1 + sqrt(1+8*DCT_fact*M)));
M_DCT = 0.5*K_DCT*(K_DCT+1);
M_rand = M - M_DCT;

fprintf('Using %d DCT and %d noiselet coefficients (%2.3f %% and %2.3f %%)\n', ...
    M_DCT, M_rand, 100*M_DCT/N, 100*M_rand/N)

vid2 = zeros(size(vid{1}, 1), size(vid{1}, 2), T_s);
vid3 = zeros(M, 1, T_s);

% Noise options
OM = zeros(M_rand, T_s);
MEAS_FUN = cell(T_s, 1);
for kk = 1:T_s
    % Set up noiselets
    q = randperm(N)';      % makes column vector of random integers 1:N
    OM(:, kk) = q(1:M_rand);    % vector of random subset of integers 1:N
    Phi  = @(z) A_noiselet (z, OM(:, kk));
    Phit = @(z) At_noiselet(z, OM(:, kk), N); 
    meas_rand.Phi = @(z) Phi(reshape(z, [], 1));
    meas_rand.Phit = @(z) reshape(Phit(z), sqrt(N), sqrt(N));

    % DCT measurements
    if M_DCT>0
        meas_dct.Phi = @(z) samp_DCT(z, M_DCT);
        meas_dct.Phit = @(z) samp_DCT_trans(z, M_DCT, [sqrt(N), sqrt(N)]);
    elseif M_DCT == 0
        meas_dct.Phi = @(z) [];
        meas_dct.Phit = @(z) 0;
    else
        error('How did you get here? Go away.')
    end
    meas_func.Phi  = @(z) [meas_dct.Phi(z); meas_rand.Phi(z)];
    meas_func.Phit = @(z) meas_dct.Phit(z(1:M_DCT)) + meas_rand.Phit(z(M_DCT+1:end));

    vid2(:, :, kk) = 2*(vid{kk}/255 - 0.5); % Renormalize to [-1, 1]
    vid3(:, :, kk) = meas_func.Phi(vid2(:, :, kk)) + sqrt(noise_var)*randn(M, 1);
    MEAS_FUN{kk} = meas_func;
end

clear vid

%if normalize_opt == 1
%    NORMS_VEC = zeros(2*N, numel(MEAS_FUN)); 
%    for LL = 1:numel(MEAS_FUN)
%        parfor PP = 1:2*N 
%            NORMS_VEC(PP, LL) = sqrt(sum(MEAS_FUN{LL}.Phi(DWTfunc.invert([zeros(PP-1, 1); 1; zeros(2*N-PP, 1)])).^2));
%            fprintf('(%d,%d);',LL,PP)
%        end
%    end
%else
%    NORMS_VEC = ones(1, T_s);
%end

varargout{1} = vid3;
varargout{2} = MEAS_FUN;
varargout{3} = vid2;

end
