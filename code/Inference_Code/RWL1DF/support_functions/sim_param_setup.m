
if strcmp(sim_type, 'yuv_vid')
    num_trials = 1;                   % Number of trials to average
    load_params.T_s = 300;            % Number of frames
    load_params.samp_factor = 0.3;    % Total sampling rate
    load_params.DCT_fact = 0;         % fraction of DCT samples
    load_params.noise_var = 0.0001;   % Measurement Noise Variance

    % Uncomment to use 4 level Daubechies Wavelet Transform
    % XFM = Wavelet('Daubechies',4,4);	
    % DWTfunc.apply = @(q) reshape(XFM*q, [], 1);
    % DWTfunc.invert = @(q) (XFM')*reshape(q, sqrt(N), sqrt(N));

    % Uncomment to use 4 level Dual-Tree Discrete Wavelet Transforms
    dwt_opts.J = 4;
    DWTfunc.apply = @(q) DTDWT2D_mat(q, dwt_opts.J);
    DWTfunc.invert = @(q) iDTDWT2D_mat(q, dwt_opts.J);

    % Set dynamics function to identity
    f_dyn = @(z) z;
    DYN_FUN{1} = f_dyn;
 
    filename = '/home/adam/Desktop/GT_Work/CodeBase_v1_0/Data/Videos/revideosequence/foreman.qcif';
    %filename = '/homeold/adam/Versioned/2011_RWL1DCS/trunk/Video_Tests/foreman.qcif';
elseif strcmp(sim_type, 'mri')

    % Choose video file (e.g. foreman in .qcif format)
    folder_name = '/home/adam/Desktop/Versioned/2011_RWL1DCS/trunk/MRIdata/';
    dataset_name = sprintf('conventional_reconstructed_images_frame%d.bmp', 1);
    filename = [folder_name,dataset_name];

    num_trials = 1;                   % Number of trials to average
    load_params.T_s = 10;             % Number of frames to analyze
    load_params.samp_factor = 0.15;   % Sampling rate
    load_params.noise_var = 0.001;   % Measurement Noise Variance

    % General Parameters
    % Uncomment to use 4 level Daubechies Wavelet Transform
    XFM = Wavelet('Daubechies',4,4);	
    DWTfunc.apply = @(q) reshape(XFM*q, [], 1);
    DWTfunc.invert = @(q) (XFM')*reshape(q, sqrt(numel(q)), sqrt(numel(q)));
    
    % Set dynamics function to identity
    f_dyn = @(z) z;
    DYN_FUN{1} = f_dyn;

elseif strcmp(sim_type, 'com_chan')
    % Some simulation parameters
    load_params.M = 80;             % Number of measurements per time step
    load_params.T_s = 5;           % Number of time steps
    load_params.N = 100;             % Max channel delay
    load_params.K = 3;		         % Number of channel reflectors
    load_params.one_step = 0;        % Choose if streaming inputs or not
    load_params.num_trials = 1;      % Number of trials to average
    load_params.noise_var = 0.0001;  % Measurement Noise Variance

    % Make sinc basis
    sinc_mat = sinc_grid(N,4,1); % chan_params.ts
    sinc_mat = sinc_mat*diag(1./sqrt(sum(sinc_mat.^2,1)));
    DWTfunc.apply = @(q) ([[sinc_mat, zeros(size(sinc_mat))];
        [zeros(size(sinc_mat)), sinc_mat]]')*q;
    DWTfunc.invert = @(q) ([[sinc_mat,zeros(size(sinc_mat))];...
        [zeros(size(sinc_mat)),sinc_mat]])*q;

    % Set dynamics function to identity
    f_dyn = @(z) z;
    DYN_FUN{1} = f_dyn;

    NORMS_VEC = ones(1, T_s);

else
    error('unknown data type')
end
