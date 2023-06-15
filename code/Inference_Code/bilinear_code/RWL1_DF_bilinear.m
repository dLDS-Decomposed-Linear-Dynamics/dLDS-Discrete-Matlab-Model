function [varargout] = RWL1_DF_bilinear(varargin)

%           %%%%%%%%%%%%%%%%%      Update Help File!!  %%%%%%%%%%%%%%%%%%%%
% [coef_dcs, recon_dcs, rMSE_dcs, PSNR_dcs] = ...
%        BPDN_DF_bilinear(MEAS_SIG, MEAS_SEL, DYN_FUN, DWTfunc, param_vals, ...
%        TRUE_VID)
%
%   The inputs are:
% 
% MEAS_SIG:   Mx1xT array of the measurements for the video frames
% MEAS_FUN:   Tx1 or 1x1 cell array of the measurement functions
% DYN_FUN:    Tx1 or 1x1 cell array of the dynamics functions
% DWTfunc:    Wavelet transform (sparsity basis)
% param_vals: struct of parameter values (has fields: lambda_val (tradeoff
%             parameter for BPDN), lambda_history (tradeoff parameter
%             between prediction and data fidelity), and tol (tolerance for
%             TFOCS solver)) 
% TRUE_VID:   Sqrt(N)xSqrt(N)xT array of the true video sequence (optional,
%             to evaluate errors)
% 
%    The outputs are:
% 
% coef_dcs:  Nx1xT array of inferred sparse coefficients
% recon_dcs: Sqrt(N)xSqrt(N)xT array of the recovered video sequence
% rMSE_dcs:  Tx1 array of rMSE values for the recovered video
% PSNR_dcs:  Tx1 array of PSNR values for the recovered video
% 
%
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated August 21, 2014. 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse Inputs
MEAS_SIG = varargin{1};
MEAS_FUN = varargin{2};
DYN_FUN = varargin{3};
DWTfunc = varargin{4};
param_vals = varargin{5};

if nargin > 5
    rMSE_calc_opt = 1;
    TRUE_VID = varargin{6};
else
    rMSE_calc_opt = 0;
end

if isfield(param_vals, 'lambda_val')
    lambda_val = param_vals.lambda_val;
else
    lambda_val = 0.001;
end
if isfield(param_vals, 'gamma_val')
    gamma_val = param_vals.gamma_val;
else
    gamma_val = 0.001;
end
if isfield(param_vals, 'alpha_a')
    alpha_a = param_vals.alpha_a;
else
    alpha_a = 0.2;
end
if isfield(param_vals, 'beta_a')
    beta_a = param_vals.beta_a;
else
    beta_a = 0.2;
end
if isfield(param_vals, 'xi_a')
    xi_a = param_vals.xi_a;
else
    xi_a = 0.9;
end
if isfield(param_vals, 'alpha_b')
    alpha_b = param_vals.alpha_b;
else
    alpha_b = 0.2;
end
if isfield(param_vals, 'beta_b')
    beta_b = param_vals.beta_b;
else
    beta_b = 0.2;
end

if isfield(param_vals, 'tol')
    TOL = param_vals.tol;
else
    TOL = 0.01;
end
if isfield(param_vals, 'num_em')
    num_em = param_vals.num_em;
else
    num_em = 5;
end
if isfield(param_vals, 'verbose')
    verbose_flag = param_vals.verbose;
else
    verbose_flag = 0;
end

DWT_apply = DWTfunc.apply;
DWT_invert = DWTfunc.invert;

meas_func = MEAS_FUN{1};
Phit = meas_func.Phit; 

M = numel(MEAS_SIG(:, :, 1));
N2 = numel(DWT_apply(Phit(MEAS_SIG(:, :, 1))));
N_b = numel(DYN_FUN);

opts.tol = TOL;
opts.printEvery = 0;
num_frames = size(MEAS_SIG, 3);
if isreal(MEAS_SIG)
    opt_set = 'R2R';
elseif ~isreal(MEAS_SIG)
    opt_set = 'C2C';
end

%% Get parameters as they correspond to the statistical model
% Using the equations:
% 
% alpha_a = (alpha+1)/lambda; xi_a = 1/(xi*lambda);     lambda_val = lambda
% beta_a = eta/(xi*lambda);   alpha_b = (tau+1)/gamma;  gamma_val = gamma
% beta_b = 1/(theta*gamma);
% 
model_alpha = alpha_a*lambda_val - 1;   
model_eta = beta_a/xi_a;
model_xi = 1/(lambda_val*xi_a);
% These parameters are not needed in their raw form
% model_theta = 1/(beta_b*gamma_val);  
% model_tau = gamma_val*alpha_val-1;

% Set parameters for estimation of the dynamics coefficients
bparams.w1 = 1/(model_alpha*model_xi);
bparams.w2 = gamma_val/model_alpha;
bparams.eta = model_eta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize coefficient and weight storage variables


coef_dcs = zeros(size(DWT_apply(Phit(MEAS_SIG(:, :, 1))), 1), size(DWT_apply(Phit(MEAS_SIG(:, :, 1))), 2),num_frames);
lambda_a = ones(size(coef_dcs));
recon_dcs = zeros(size(Phit(MEAS_SIG(:, :, 1)), 1), size(Phit(MEAS_SIG(:, :, 1)), 2),num_frames);
bcoef_dcs = zeros(N_b,1,num_frames);
lambda_b = ones(size(bcoef_dcs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for initial frame

meas_func = MEAS_FUN{1};
Phi  = meas_func.Phi;
Phit = meas_func.Phit; 

tic
for em_iter = 1:num_em
    % Generate A/A^t
    Af = @(x) Phi(DWT_invert(x));
    Ab = @(x) DWT_apply(Phit(x));
    A = linop_handles([M, N2], Af, Ab, opt_set);
    % Solve coefficients given weights
    coef_dcs(:, :, 1) = solver_L1RLS(A, MEAS_SIG(:, :, 1), lambda_val, zeros(N2, 1), opts );
    % Update weights
    lambda_a(:, :, 1) = alpha_a./(abs(coef_dcs(:, :, 1)) +  beta_a);
    
end
TIME_ITER = toc; % Get timing

% Save/display reconstruction results
recon_dcs(:, :, 1) = DWT_invert(real(coef_dcs(:, :, 1)));
if rMSE_calc_opt == 1
    rMSE_dcs(1) = sum(sum((recon_dcs(:, :, 1) - TRUE_VID(:, :, 1)).^2))/sum(sum(TRUE_VID(:, :, 1).^2));
    PSNR_dcs(1) = psnr(real(recon_dcs(:, :, 1)), TRUE_VID(:, :, 1), 2);
    if verbose_flag == 1
        fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', 1, num_frames, TIME_ITER, PSNR_dcs(1), rMSE_dcs(1))
    end
else
    if verbose_flag == 1
        fprintf('Finished frame %d of %d in %f seconds. \n', 1, num_frames, TIME_ITER)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for rest of frames
for kk = 2:num_frames
    tic
    % Four Steps:
    % M-step 1: find signal coefficients given first variance parameters
    % M-step 2: fing dynamics coefficients given second variance parameters
    % E-step 1: calculate new first variance parameters given signal and dynamics coefficients
    % E-step 2: calcualte new second variance parameters given dynamics coefficinets
	
    %% parameters that are good to calculate

    % Get the Measurement function    
    Phi  = meas_func.Phi;
    Phit = meas_func.Phit;
 
    if numel(MEAS_FUN) == 1
        meas_func = MEAS_FUN{1};
    elseif numel(MEAS_FUN) == num_frames
        meas_func = MEAS_FUN{kk};
    else
        error('You need either the same measurement function for all time or one dynamics function per time-step!')
    end

    % Project previous estimate throught he dynamics functions. The inner product of 
    % each row with the vector of dynamics coefficients gives the estimate for that 
    % element of the signal coefficients
    f_dyn = zeros(size(coef_dcs, 1),numel(DYN_FUN));
    for ll = 1:N_b
        f_dyn(:, ll) = DYN_FUN{ll}*coef_dcs(:,:,kk-1);
    end
    
    %% Start EM iterations
    for em_iter = 1:num_em
        
        %% M-step 1: Get sparse coefficients
        Af = @(x) Phi(DWT_invert(x./lambda_a(:,:,kk)));     % Set A/At for BPDN substep
        Ab = @(x) DWT_apply(Phit(x))./lambda_a(:,:,kk);
	    A = linop_handles([M, N2], Af, Ab, opt_set);
        
        coef_dcs(:, :, kk) = solver_L1RLS(A, MEAS_SIG(:, :, kk), lambda_val, zeros(N2, 1), opts);   % Optimize the BPDN objective function with TFOCS

        %% M-step 2 (the hard one - hybrid sampling/convex optimization)
        bcoef_dcs(:, :, kk) = opt_regularized_AC(f_dyn, lambda_a(:, :, kk), lambda_b(:, :, kk), bparams); 
        
        %% E-step 1: Update coefficient weights
        lambda_a(:, :, kk) = alpha_a./(abs(coef_dcs(:, :, kk)) + xi_a*abs(f_dyn*bcoef_dcs(:, :, kk)) +  beta_a); 
        
        %% E-step 2: Update dynamics weights
        lambda_b(:, :, kk) = alpha_b./(abs(bcoef_dcs(:, :, kk)) + beta_b);   
%         if (sum(isnan(lambda_b(:, :, kk))) > 0)||(sum(isnan(lambda_a(:, :, kk))) > 0)
%             sum(isnan(lambda_b(:, :, kk)))
%         end
    end
    TIME_ITER = toc;
    
    %% End EM iterations: Save/display reconstruction results
    
    recon_dcs(:, :, kk) = DWT_invert(real(coef_dcs(:, :, kk)));
    if rMSE_calc_opt == 1
        rMSE_dcs(kk) = sum(sum((recon_dcs(:, :, kk) - TRUE_VID(:, :, kk)).^2))/sum(sum(TRUE_VID(:, :, kk).^2));
        PSNR_dcs(kk) = psnr(real(recon_dcs(:, :, kk)), TRUE_VID(:, :, kk), 2);
        if verbose_flag == 1
            fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', kk, num_frames, TIME_ITER, PSNR_dcs(kk), rMSE_dcs(kk))
        end
    else
        if verbose_flag == 1
            fprintf('Finished frame %d of %d in %f seconds.\n', kk, num_frames, TIME_ITER)
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set ouptputs

if (rMSE_calc_opt == 1)
    if nargout > 0
        varargout{1} = {coef_dcs, bcoef_dcs, lambda_a, lambda_b};
    end
    if nargout > 1
        varargout{2} = recon_dcs;
    end
    if nargout > 2
        varargout{3} = rMSE_dcs;
    end
    if nargout > 3
        varargout{4} = PSNR_dcs;
    end
    if nargout > 4
        for kk = 5:nargout
            varargout{kk} = [];
        end
    end
elseif (rMSE_calc_opt ~= 1)
    if nargout > 0
        varargout{1} = {coef_dcs, bcoef_dcs, lambda_a, lambda_b};
    end
    if nargout > 1
        varargout{2} = recon_dcs;
    end
    if nargout > 2
        for kk = 3:nargout
            varargout{kk} = [];
        end
    end
else
    error('How did you get here?')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
