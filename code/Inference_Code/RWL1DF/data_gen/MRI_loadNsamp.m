function varargout = MRI_loadNsamp(varargin)


file_name = varargin{1};
param_obj = varargin{2};


if (~isfield(param_obj,'samp_factor'))||(~isfield(param_obj,'T_s'))
    error('Need to specify samp_factor and T_s')
end

num_frames = param_obj.T_s;
samp_factor = param_obj.samp_factor;

if ~isfield(param_obj,'samp_opt')
    samp_opt = 'Cols';
else
    samp_opt = param_obj.samp_opt;
    if isempty(DCT_fact)
        samp_opt = 'Cols';
    end
end

if ~isfield(param_obj,'n_var')
    n_var = 0.0001;
else
    n_var = param_obj.n_var;
    if isempty(n_var)
        n_var = 0.0001;
    end
end
if ~isfield(param_obj,'same_samp')
    same_samp = 0;
else
    same_samp = param_obj.same_samp;
    if isempty(same_samp)
        same_samp = 0;
    end
end

% samp_opt = 'Cols';
% n_var = 0.0001; % Variance of added noise
plot_opt = 0;
MEAS_FUN = cell(1, num_frames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data

% Set image filenames
tmp = imread(file_name);
N = size(tmp);

% Load images and zero-pad
aheart_image = zeros(N(1)+64, N(2), num_frames);
aheart_image(:, :, 1) = [zeros(32, 256); tmp; zeros(32, 256)];
clear tmp
for kk = 2:num_frames
    file_name = sprintf('conventional_reconstructed_images_frame%d.bmp', kk);
    aheart_image(:, :, kk) = [zeros(32, 256); imread(file_name); zeros(32, 256)];
end

% Re-scale to [-1,1]
aheart_image = 2*aheart_image/255 - 1;
% Image size:
N = [256, 256];

if plot_opt == 1
    % Display image sequence
    for kk = 1:num_frames
        imagesc(aheart_image(:, :, mod(kk, num_frames)+1))
        colormap gray
        axis image
        axis off
        pause(0.05)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose sampling pattern

% Make sampling matrix
if strcmp(samp_opt, 'Pxls')
    % Random pixels
    % Base Sampling
    num_samps = round(samp_factor*numel(aheart_image(:, :, 1)));
    [X,Y] = meshgrid(linspace(-1, 1, N(2)), linspace(-1, 1, N(1)));
    pdf = exp(-(X.^2 + Y.^2)/(0.4));
    pdf = num_samps*pdf/sum(sum(pdf));
    iter = 40;
    tol = 50;
    [mask_mult,stat,actpctg] = genSampling(pdf,iter,tol);
    % Per Image Sampling
    mask_mult_cell = cell(num_frames, 1);
    for ll = 1:num_frames
	if same_samp == 1
		mask_mult_cell{ll} = mask_mult;
	else
		[mask_mult_cell{ll},~,~] = genSampling(pdf,iter,tol);
	end
    end
elseif strcmp(samp_opt, 'Cols')
    % Random Columns
    % Base Sampling
    num_samps = round(samp_factor*numel(aheart_image(:, 1, 1)));
    X = linspace(-1, 1, N(1));
    pdf = exp(-abs(X)/(0.2));
    pdf = num_samps*pdf/sum(sum(pdf));
    iter = 50;
    tol = 1;
    [mask_mult,stat,actpctg] = genSampling(pdf,iter,tol);
    mask_mult = ones(256, 1)*mask_mult;
    % Per Image Sampling
    mask_mult_cell = cell(num_frames, 1);
    for ll = 1:num_frames
	if same_samp == 1
		mask_mult_cell{ll} = mask_mult;
	else
		[mask_mult_cell{ll},~,~] = genSampling(pdf,iter,tol);
		mask_mult_cell{ll} = ones(256, 1)*mask_mult_cell{ll};
	end
    end
elseif strcmp(samp_opt, 'Brain')
    % From Brain image
    load brain512
    % Base Sampling
    mask_mult = mask(1:2:end, 1:2:end);
    mask_mult_cell = cell(num_frames, 1);
    for ll = 1:num_frames
	mask_mult_cell{ll} = mask(1:2:end, 1:2:end);
    end
else
    error('Unknown Sampling Scheme!')
end

if plot_opt == 1
    % Display sampling pattern
    imagesc(mask_mult)
    colormap gray
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take FFTs of data

% Pre-allocate arrays to hold FFT data
% aheart_fft = zeros(size(aheart_image(:, :, kk)));
% aheart_fft_sub = zeros(size(aheart_image(:, :, kk)));
aheart_fft_sub2 = zeros(numel(aheart_image(:, :, 1)), 1, num_frames);

% FT = p2DFT(mask_mult, N, 1, 2);

for kk = 1:num_frames
    FT2 = p2DFT(mask_mult_cell{kk}, size(aheart_image(:, :, 1)), 1, 2);
    meas_func.Phi  = @(z) reshape(FT2*z, [], 1); 
    meas_func.Phit = @(z) (FT2')*reshape(z,size(aheart_image(:, :, 1)));
    MEAS_FUN{kk} = meas_func;

    % Take partial FFTs and add noise
%    aheart_fft_sub(:, :, kk) = FT*aheart_image(:, :, kk);
%    aheart_fft_sub(:, :, kk) = aheart_fft_sub(:, :, kk) + ...
%	sqrt(n_var/2)*(randn(size(aheart_fft_sub(:, :, kk))) + ...
%	1j*randn(size(aheart_fft_sub(:, :, kk))));
    aheart_fft_sub2(:, :, kk) = meas_func.Phi(aheart_image(:, :, kk));
    aheart_fft_sub2(:, :, kk) = aheart_fft_sub2(:, :, kk) + ...
	sqrt(n_var/2)*(randn(numel(aheart_fft_sub2(:, :, kk)),1) + ...
	1j*randn(numel(aheart_fft_sub2(:, :, kk)),1));
end

%if plot_opt == 1
%    % Plot subsampled FFTs
%    for kk = 1:num_frames
%	subplot(1, 2, 1), imagesc(20*log10(abs(aheart_fft_sub(:, :, kk))))
%	colormap gray
%	axis image
%	subplot(1, 2, 2), imagesc(20*log10(abs(aheart_fft_sub2(:, :, kk))))
%	colormap gray
%	axis image
%	pause(0.05)
%    end
%end

varargout{1} = aheart_fft_sub2;
varargout{2} = MEAS_FUN;
varargout{3} = aheart_image;

end
