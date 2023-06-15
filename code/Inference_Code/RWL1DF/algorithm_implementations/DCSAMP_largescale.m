function [varargout] = DCSAMP_video(varargin)

% [vid_coef_amp, vid_recon_amp, vid_rMSE_amp, vid_PSNR_amp] = ...
%        DCSAMP_video(MEAS_SIG, MEAS_SEL, DWTfunc, param_vals, ...
%        TRUE_VID)
%
% This function is essentially a wrapper for Justin Ziniel and Phil
% Schniter's DCS-AMP filtering algorithm. 
% 
%   The inputs are:
% 
% MEAS_SIG:   Mx1xT array of the measurements for the video frames
% MEAS_FUN:   Tx1 or 1x1 cell array of the measurement functions
% DWTfunc:    Wavelet transform (sparsity basis)
% param_vals: struct of parameter values (has fields: lambda_val (tradeoff
%             parameter for BPDN), lambda_history (tradeoff parameter
%             between prediction and data fidelity), and tol (tolerance for
%             TFOCS solver)) 
% NORMS_VEC:  NxT array In the case where an implicit measurement matrix is
%             not normalized, this matrix should contain the norms of each
%             column of the measurement matrix (Phi(DWT_inverse(..))) in
%             order to manually normalize the columns. If left empty, this
%             will default to ones.
% DCS_options: Struct containing options for DCS_AMP tp run (i.e. number of
%             smoothing iterations, choice of updating, etc)
% TRUE_VID:   Sqrt(N)xSqrt(N)xT array of the true video sequence (optional,
%             to evaluate errors)
% 
%    The outputs are:
% 
% vid_coef_amp:  Nx1xT array of inferred sparse coefficients
% vid_recon_amp: Sqrt(N)xSqrt(N)xT array of the recovered video sequence
% vid_rMSE_amp:  Tx1 array of rMSE values for the recovered video
% vid_PSNR_amp:  Tx1 array of PSNR values for the recovered video
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

MEAS_SIG = varargin{1};
MEAS_FUN = varargin{2};
DWTfunc = varargin{3};
DCS_params = varargin{4};
NORMS_VEC = varargin{5};
if isempty(NORMS_VEC)
    NORMS_VEC = ones(1, numel(MEAS_FUN));
end

if nargin > 5
    DCS_options = varargin{6};
else
    DCS_options.smooth_iter = -1;
    DCS_options.update = 1;
end


if nargin > 6
    rMSE_calc_opt = 1;
    TRUE_VID = varargin{7};
else
    rMSE_calc_opt = 0;
end

if nargin > 7
    rMSE_plot_opt = varargin{7};
else
    rMSE_plot_opt = 0;
end

meas_func = MEAS_FUN{1};
Phit = meas_func.Phit;
DWT_apply = DWTfunc.apply;

N2 = numel(DWT_apply(Phit(MEAS_SIG(:, :, 1))));
DCS_params.lambda_0 = DCS_params.lambda_0*ones(N2, 1);
DCS_params.eta_0    = DCS_params.eta_0*ones(N2, 1);    	
DCS_params.kappa_0  = DCS_params.kappa_0*ones(N2, 1);
DCS_params.alpha    = DCS_params.alpha*ones(N2, 1);
DCS_params.rho      = DCS_params.rho*ones(N2, 1);


MEAS_FUN_DCS = cell(size(MEAS_SIG, 3), 1);
y_DCS = cell(size(MEAS_SIG, 3), 1);

if numel(MEAS_FUN) == 1
    for kk = 1:size(MEAS_SIG, 3)
        MEAS_FUN_DCS{kk} = @(z, m) DCS_meas_func(z, m, MEAS_FUN{1}, DWTfunc, NORMS_VEC);
        y_DCS{kk} = MEAS_SIG(:, :, kk);
    end
else
    for kk = 1:size(MEAS_SIG, 3)
        MEAS_FUN_DCS{kk} = @(z, m) DCS_meas_func(z, m, MEAS_FUN{kk}, DWTfunc, NORMS_VEC(:, kk));
        y_DCS{kk} = MEAS_SIG(:, :, kk);
    end
end

for kk = 1:size(MEAS_SIG, 3)
    y_DCS{kk} = MEAS_SIG(:, :, kk);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Algorithm and Collect Results

tic
[x_hat, ~, ~] = sp_multi_frame_fxn(y_DCS, MEAS_FUN_DCS, DCS_params, DCS_options);
toc

if rMSE_calc_opt == 1
    for kk = 1:size(TRUE_VID, 3)
        vid_coef_amp(:, :, kk) = x_hat{kk};
        vid_recon_amp(:, :, kk) = DWTfunc.invert((1./NORMS_VEC(:, kk)).*vid_coef_amp(:, :, kk));
        vid_rMSE_amp(kk) = sum(sum((vid_recon_amp(:, :, kk) - TRUE_VID(:, :, kk)).^2))/sum(sum(TRUE_VID(:, :, kk).^2));
        vid_PSNR_amp(kk) = psnr(real(vid_recon_amp(:, :, kk)), TRUE_VID(:, :, kk), 2);
        
        if rMSE_plot_opt == 1
            plot(vid_rMSE_amp, 'LineWidth', 3)
            xlabel('Frame Number', 'FontSize', 22)
            ylabel('rMSE', 'FontSize', 22)
            set(gca, 'FontSize', 20)

            disp(mean(vid_rMSE_amp))
            disp(median(vid_rMSE_amp))
        else
            % Do NOTHING
        end
    end
else
    for kk = 1:size(TRUE_VID, 3)
        vid_coef_amp(:, :, kk) = x_hat{kk};
        vid_recon_amp(:, :, kk) = DWTfunc.invert((1./NORMS_VEC(:, kk)).*vid_coef_amp(:, :, kk));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set ouptputs

if (rMSE_calc_opt == 1)
    if nargout > 0
        varargout{1} = vid_coef_amp;
    end
    if nargout > 1
        varargout{2} = vid_recon_amp;
    end
    if nargout > 2
        varargout{3} = vid_rMSE_amp;
    end
    if nargout > 3
        varargout{4} = vid_PSNR_amp;
    end
    if nargout > 4
        for kk = 5:nargout
            varargout{kk} = [];
        end
    end
elseif (rMSE_calc_opt ~= 1)
    if nargout > 0
        varargout{1} = vid_coef_amp;
    end
    if nargout > 1
        varargout{2} = vid_recon_amp;
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
