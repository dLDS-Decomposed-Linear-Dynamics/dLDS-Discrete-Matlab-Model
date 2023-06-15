function [varargout] = RWL2mod_largescale(varargin)

% [vid_coef_rwcs, vid_recon_rwcs, vid_rMSE_rwcs, vid_PSNR_rwcs] = ...
%        RWL2mod_largescale(MEAS_SIG, MEAS_FUN, DYN_FUN, DWTfunc, ...
%        param_vals, TRUE_VID)
%
%   The inputs are:
% 
% MEAS_SIG:   Mx1xT array of the measurements for the video frames
% MEAS_FUN:   Tx1 or 1x1 cell array of the measurement functions
% DYN_FUN:    Tx1 or 1x1 cell array of the dynamics functions
% DWTfunc:    Wavelet transform (sparsity basis)
% param_vals: struct of parameter values (has fields: lambda_val (tradeoff
%             parameter for BPDN), rwl1_mult (multiplicative parameter for
%             the re-weighting), rwl1_reg (additive value in denominator of
%             the re-weighting), beta (multiplicative value in denominator 
%             of the re-weighting) and tol (tolerance for TFOCS solver))
% TRUE_VID:   Sqrt(N)xSqrt(N)xT array of the true video sequence (optional,
%             to evaluate errors)
% 
%    The outputs are:
% 
% vid_coef_rwcs:  Nx1xT array of inferred sparse coefficients
% vid_recon_rwcs: Sqrt(N)xSqrt(N)xT array of the recovered video sequence
% vid_rMSE_rwcs:  Tx1 array of rMSE values for the recovered video
% vid_PSNR_rwcs:  Tx1 array of PSNR values for the recovered video
% 
% 
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated August 20, 2012. 
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

if nargin > 6
    NORM_vals = varargin{7};
else
    NORM_vals = ones(1, numel(MEAS_FUN));
end

DWT_apply = DWTfunc.apply;
DWT_invert = DWTfunc.invert;

if isfield(param_vals, 'lambda_val')
    lambda_val = param_vals.lambda_val;
else
    lambda_val = 0.001;
end
if isfield(param_vals, 'rwl1_reg')
    rwl1_reg = param_vals.rwl1_reg;
else
    rwl1_reg = 0.2;
end
if isfield(param_vals, 'rwl1_mult')
    rwl1_mult = param_vals.rwl1_mult;
else
    rwl1_mult = 0.4;
end
if isfield(param_vals, 'beta')
    w1 = param_vals.beta;
else
    w1 = 1;
end
if isfield(param_vals, 'tol')
    TOL = param_vals.tol;
else
    TOL = 0.001;
end

meas_func = MEAS_FUN{1};
Phit = meas_func.Phit; 

M = numel(MEAS_SIG(:, :, 1));
N2 = numel(DWT_apply(Phit(MEAS_SIG(:, :, 1))));

num_frames = size(MEAS_SIG, 3);
opts.tol = TOL;
opts.maxIts = 1000;
opts.printEvery = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run RWL1 on each frame

% Run first time-step
kk = 1;
tic
meas_func = MEAS_FUN{1};
Phi  = meas_func.Phi;
Phit = meas_func.Phit; 

weights = 1;
for nn = 1:10
    % M Step: Solve BPDN
    Af_tmp = @(x) Phi(DWT_invert((1./NORM_vals(:, 1)).*(x)));
    Ab_tmp = @(x) (1./NORM_vals(:, 1)).*DWT_apply(Phit(x));
    
    Af = @(x) Af_tmp(Ab_tmp(x)./weights);
    Ab = @(x) Af(x);
    A  = linop_handles([M, M], Af, Ab, 'R2R');
    
    res = tfocs(smooth_quad, { A, -reshape(MEAS_SIG(:, :, 1), [], 1) }, [], [], opts);
    res = Ab_tmp(res)./weights;
    
    % E-step: Update weights
    weights = 1./(abs(real(res)).^2 + rwl1_reg);

    if rMSE_calc_opt == 1
        im_res = DWT_invert(res./NORM_vals(:, 1));
        temp_rMSE = sum(sum((im_res - TRUE_VID(:, :, 1)).^2))/sum(sum(TRUE_VID(:, :, 1).^2));
        fprintf('Finished RW iteration %d, rMSE = %f.. \n', nn, temp_rMSE)
    else
        fprintf('Finished RW iteration %d.\n', nn)
    end
end

if rMSE_calc_opt == 1
    % Save results
    vid_coef_drl2(:, :, 1) = res;
    vid_recon_drl2(:, :, 1) = im_res;
    vid_rMSE_drl2(1) = sum(sum((vid_recon_drl2(:, :, 1) - TRUE_VID(:, :, 1)).^2))/sum(sum(TRUE_VID(:, :, 1).^2));
    vid_PSNR_drl2(1) = psnr(vid_recon_drl2(:, :, 1), TRUE_VID(:, :, 1), 2);
    TIME_ITER = toc;
    fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', 1, num_frames, TIME_ITER, vid_PSNR_drl2(1), vid_rMSE_drl2(1))
else
    vid_coef_drl2(:, :, 1) = res;
    vid_recon_drl2(:, :, 1) = im_res;
    TIME_ITER = toc;
    fprintf('Finished frame %d of %d in %f seconds.\n', 1, num_frames, TIME_ITER)
end

% Solve the remaining frames
for kk = 2:num_frames
    tic
    if numel(DYN_FUN) == 1
        f_dyn = DYN_FUN{1};
    elseif numel(DYN_FUN) == num_frames
        f_dyn = DYN_FUN{kk};
    else
        error('You need either the same dynamics function for all time or one dynamics function per time-step!')
    end
    if numel(MEAS_FUN) == 1
        meas_func = MEAS_FUN{1};
    elseif numel(MEAS_FUN) == num_frames
        meas_func = MEAS_FUN{kk};
    else
        error('You need either the same dynamics function for all time or one dynamics function per time-step!')
    end
    
    % Predict weights
    x_pred = NORM_vals(:, kk).*DWT_apply(f_dyn(vid_recon_drl2(:, :, kk-1)));
    % Set measurement function for the time-step
    Phi  = meas_func.Phi;
    Phit = meas_func.Phit;
    % Initialize EM weights
    weights = ones(size(x_pred));
    for nn = 1:10
        % M-Step: Solve BPDN
        
        Af_tmp = @(x) Phi(DWT_invert((1./NORM_vals(:, kk)).*x));
        Ab_tmp = @(x) (1./NORM_vals(:, kk)).*DWT_apply(Phit(x));
        Af = @(x) Af_tmp(Ab_tmp(x)./weights);
        Ab = @(x) Af(x);
        A  = linop_handles([M, M], Af, Ab, 'R2R');
        
        res = tfocs(smooth_quad, {A, -reshape(MEAS_SIG(:, :, 1), [], 1)}, [], [], opts);
        res = Ab_tmp(res)./weights;
        
        % E-Step: Update weights
	supp_est = abs(vid_coef_drl2(:, :, kk-1)) >= 0.01;
	%supp_est = abs(x_pred) >= 0.01;
        weights(supp_est) = 1./(w1*(abs(real(res(supp_est))).^2) + rwl1_reg);
        weights(~supp_est) = 1./((abs(real(res(~supp_est))).^2) + rwl1_reg);
        
        if rMSE_calc_opt == 1
            im_res = DWT_invert(res./NORM_vals(:, kk));
            temp_rMSE = sum(sum((im_res - TRUE_VID(:, :, kk)).^2))/sum(sum(TRUE_VID(:, :, kk).^2));
            fprintf('Finished RW iteration %d, rMSE = %f, support = %f.. \n', nn, temp_rMSE, sum(supp_est)/numel(supp_est))
        else
            fprintf('Finished RW iteration %d.\n', nn)
        end
    end
    
    if rMSE_calc_opt == 1
        vid_coef_drl2(:, :, kk) = res;
        vid_recon_drl2(:, :, kk) = im_res;
        vid_rMSE_drl2(kk) = sum(sum((vid_recon_drl2(:, :, kk) - TRUE_VID(:, :, kk)).^2))/sum(sum(TRUE_VID(:, :, kk).^2));
        vid_PSNR_drl2(kk) = psnr(real(vid_recon_drl2(:, :, kk)), TRUE_VID(:, :, kk), 2);
        TIME_ITER = toc;
        fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', kk, num_frames, TIME_ITER, vid_PSNR_drl2(kk), vid_rMSE_drl2(kk))
    else
        im_res = DWT_invert(res./NORM_vals(:, kk));
        vid_coef_drl2(:, :, kk) = res;
        vid_recon_drl2(:, :, kk) = im_res;
        TIME_ITER = toc;
        fprintf('Finished frame %d of %d in %f seconds.\n', kk, num_frames, TIME_ITER)
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set ouptputs

if (rMSE_calc_opt == 1)
    if nargout > 0
        varargout{1} = vid_coef_drl2;
    end
    if nargout > 1
        varargout{2} = vid_recon_drl2;
    end
    if nargout > 2
        varargout{3} = vid_rMSE_drl2;
    end
    if nargout > 3
        varargout{4} = vid_PSNR_drl2;
    end
    if nargout > 4
        for kk = 5:nargout
            varargout{kk} = [];
        end
    end
else
    if nargout > 0
        varargout{1} = vid_coef_drl2;
    end
    if nargout > 1
        varargout{2} = vid_recon_drl2;
    end
    if nargout > 2
        for kk = 3:nargout
            varargout{kk} = [];
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
