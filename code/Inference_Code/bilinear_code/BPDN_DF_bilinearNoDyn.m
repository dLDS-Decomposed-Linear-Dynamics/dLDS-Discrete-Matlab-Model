function [varargout] = BPDN_DF_bilinearNoDyn(varargin)

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
% Last updated August 21, 2012. 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse Inputs

MEAS_SIG   = varargin{1};
MEAS_FUN   = varargin{2};
DYN_FUN    = varargin{3};
DWTfunc    = varargin{4};
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
if isfield(param_vals, 'lambda_history')
    lambda_history = param_vals.lambda_history;
else
    lambda_history = 0.2;
end
if isfield(param_vals, 'lambda_b')
    lambda_b = param_vals.lambda_b;
else
    lambda_b = 0.2;
end
if isfield(param_vals, 'lambda_historyb')
    lambda_historyb = param_vals.lambda_historyb;
else
    lambda_historyb = 0;
end
if isfield(param_vals, 'tol')
    TOL = param_vals.tol;
else
    TOL = 0.01;
end
if isfield(param_vals, 'verbose')
    verbose_flag = param_vals.verbose;
else
    verbose_flag = 0;
end



DWT_apply  = DWTfunc.apply;
DWT_invert = DWTfunc.invert;

meas_func = MEAS_FUN{1};
Phit      = meas_func.Phit; 

M   = numel(MEAS_SIG(:, :, 1));
N2  = numel(DWT_apply(Phit(MEAS_SIG(:, :, 1))));
N_b = numel(DYN_FUN);

if isfield(param_vals,'x_st')
    x_st = param_vals.x_st;
else
    x_st = zeros(N2,1);
end

opts.tol        = TOL;
opts.printEvery = 0;
num_frames      = size(MEAS_SIG, 3);
if isreal(MEAS_SIG)
    opt_set = 'R2R';
elseif ~isreal(MEAS_SIG)
    opt_set = 'C2C';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for initial frame

meas_func = MEAS_FUN{1};
Phi       = meas_func.Phi;
Phit      = meas_func.Phit; 

Af = @(x) Phi(DWT_invert(x));
Ab = @(x) DWT_apply(Phit(x));
A  = linop_handles([M, N2], Af, Ab, opt_set);

bcoef_dcs(:, :, 1) = zeros(N_b,1);

for kk = 1:num_frames
    res    = solver_L1RLS(A, MEAS_SIG(:, :, kk), lambda_val, x_st, opts );
    im_res = DWT_invert(real(res));

    % Save reconstruction results
    coef_dcs(:, :, kk)  = res;
    recon_dcs(:, :, kk) = im_res;
    tic
    if rMSE_calc_opt == 1
        rMSE_dcs(1) = sum(sum((recon_dcs(:, :, 1) - TRUE_VID(:, :, 1)).^2))/sum(sum(TRUE_VID(:, :, 1).^2));
        PSNR_dcs(1) = psnr(real(recon_dcs(:, :, 1)), TRUE_VID(:, :, 1), 2);
        TIME_ITER   = toc;
        if verbose_flag == 1
            fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', 1, num_frames, TIME_ITER, PSNR_dcs(1), rMSE_dcs(1))
        end
    else
        TIME_ITER = toc;
        if verbose_flag == 1
            fprintf('Finished frame %d of %d in %f seconds. \n', 1, num_frames, TIME_ITER)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set ouptputs

if (rMSE_calc_opt == 1)
    if nargout > 0
        varargout{1} = {coef_dcs, bcoef_dcs};
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
        varargout{1} = {coef_dcs, bcoef_dcs};
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
