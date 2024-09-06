function [varargout] = BPDN_DF_bilinear_behavior(varargin)

% [coef_dcs, recon_dcs, rMSE_dcs, PSNR_dcs] = ...
%        BPDN_DF_bilinear(MEAS_SIG, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, ...
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

Psi        = varargin{6};
behavior   = varargin{7};

% disp(size(MEAS_SIG))
% disp(size(MEAS_FUN{1}))
% disp(size(Psi))
% disp(size(behavior))

if nargin > 5
    ifBehavior = 1;
end

if nargin > 7
    rMSE_calc_opt = 1;
    TRUE_VID = varargin{8};
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

if isfield(param_vals, 'solver_type')
    if strcmp(param_vals.solver_type, 'cvx')
        doCVX = true;
        doTFOCS = false;
    elseif strcmp(param_vals.solver_type, 'tfocs')
        doCVX = false;
        doTFOCS = true;
    else
        % defaults to fista
        doCVX = false;
        doTFOCS = false; 
    end
else
    % defaults to fista
    doCVX = false;
    doTFOCS = false; 
end

if isfield(param_vals, 'CVX_Precision')
    CVX_Precision = param_vals.CVX_Precision;
else
    CVX_Precision = 'default';
end

if isfield(param_vals, 'deltaDynamics');    deltaOpt = param_vals.deltaDynamics;
else;                                       deltaOpt = false;
end

if isfield(param_vals, 'debias');    debiasOpt = param_vals.debias;
else;                                debiasOpt = true; % debias by default
end

if isfield(param_vals, 'lambda_behavior')
    lambda_behavior = param_vals.lambda_behavior;
else
    lambda_behavior = 0.1;
end

if isfield(param_vals, 'nBhv')
    nBhv = param_vals.nBhv;
else
    nBhv  = size(Psi,2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some initializations and setups

DWT_apply  = DWTfunc.apply;
DWT_invert = DWTfunc.invert;

meas_func = MEAS_FUN{1};
Phit      = meas_func.Phit; 

M   = numel(MEAS_SIG(:, :, 1));
N2  = numel(DWT_apply(Phit(MEAS_SIG(:, :, 1))));
N_b = numel(DYN_FUN);

opts.tol        = TOL;
opts.printEvery = 0;
num_frames      = size(MEAS_SIG, 3);
% disp(num_frames)
if isreal(MEAS_SIG);       opt_set = 'R2R';   % Set real vals if needed
elseif ~isreal(MEAS_SIG);  opt_set = 'C2C';   % Set imag vals if needed
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for initial frame

meas_func = MEAS_FUN{1};
Phi       = meas_func.Phi;
Phit      = meas_func.Phit; 

Af = @(x) Phi(DWT_invert(x));
Ab = @(x) DWT_apply(Phit(x));
A  = linop_handles([M, N2], Af, Ab, opt_set);

Psi       = meas_func.Psi;
Psit      = meas_func.Psit; 

% disp(size(Psi))

if ifBehavior
    APsif = @(x) Psi(DWT_invert(x));
    APsib = @(x) DWT_apply(Psit(x));
end

res    = solver_L1RLS(A, MEAS_SIG(:, :, 1), lambda_val, zeros(N2, 1), opts );
im_res = DWT_invert(real(res));

% Save reconstruction results
coef_dcs(:, :, 1)  = res;
bcoef_dcs(:,:,1)   = zeros(N_b,1);
recon_dcs(:, :, 1) = im_res;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for rest of frames

% whatIsD = DWT_apply(Phit(eye(N2))).';
whatIsD = Af(eye(N2));
if ifBehavior
   whatIsPsi = APsif(eye(N_b));
end


for kk = 2:num_frames
    tic
    
    % Get the Dynamics function
    f_dyn = zeros(size(coef_dcs, 1),numel(DYN_FUN));
    for ll = 1:N_b
        f_dyn(:, ll) = DYN_FUN{ll}*coef_dcs(:,:,kk-1);
    end
    
    % Get the Measurement function
    if numel(MEAS_FUN) == 1;                  meas_func = MEAS_FUN{1};     % If one meas fun, use the same everywhere
    elseif numel(MEAS_FUN) == num_frames;     meas_func = MEAS_FUN{kk};    % Otherwise update at each time-point
    else;  error('You need either the same measurement function for all time or one dynamics function per time-step!')
    end

    % Set  up A and At for TFOCS
    Phi  = meas_func.Phi;
    Phit = meas_func.Phit;
        
    
    %fista
    if deltaOpt;  yNow = MEAS_SIG(:,:,kk) - MEAS_SIG(:,:,kk-1);
    else;         yNow = MEAS_SIG(:,:,kk);                        end
        
    if ifBehavior
        
        if lambda_historyb == 0
            % disp('lambda_historyb == 0')
            % disp(size(yNow))
            % disp(size(whatIsD))
            % disp(size(whatIsPsi))
            % disp(N_b)
            % disp(N2)
            % disp(bcoef_dcs(:,:,kk-1))
            % disp(size(zeros(size(whatIsPsi,1),N2)))
            % disp(size(whatIsPsi))
            % disp(size(behavior))
            % disp('how to invert psi?')
            % disp(size(whatIsD)) % size: data dim, latent dim
            % disp(size(zeros(size(whatIsD,1),N_b))) % size: data dim, nF
	    % disp('FIXME: something is still wrong with Psi and behavior')
            % disp(size(zeros(size(Psi,1),N2))) % size: nBhv, latent dim
            % disp(size(bcoef_dcs(:,:)))
            % disp(size(Psi*behavior(kk,:).')) % size: nF, just this time point
            
            thisBehaviorSample = behavior{:,:}; 
            % disp(thisBehaviorSample)
            % disp(size(thisBehaviorSample(kk-1)))
            % disp(thisBehaviorSample(:,kk-1))



            A = [zeros(N_b,N2) eye(N_b,N_b); sqrt(lambda_history)*eye(N2,N2) sqrt(lambda_history)*-1*f_dyn;whatIsD zeros(size(whatIsD,1),N_b); zeros(size(whatIsPsi,1),N2) sqrt(lambda_behavior)*whatIsPsi];
            % disp(class(bcoef_dcs(:,:,kk-1)))
            % disp(class(zeros(N2,1)))
            % disp(class(yNow))
            % disp(class(thisBehaviorSample(kk-1)))
            Y = [bcoef_dcs(:,:,kk-1); zeros(N2,1); yNow; thisBehaviorSample(:,kk-1)]; % add on nBhv by 1

            lamVec = [lambda_val*ones(N2,1);lambda_b*ones(N_b,1)];

            opts.pos        = false;
            opts.lambda     = lamVec(:);
            opts.check_grad = 0;
    
            res = fista_lasso(Y, A, [], opts); %FIXME: augmented Y needed, not just measurement
       
        else
            
%             A      = [sqrt(lambda_historyb)*eye(N_b,N_b); f_dyn];
%             Y      = [sqrt(lambda_historyb)*bcoef_dcs(:,:,kk-1); yNow];
%             lamVec = lambda_b*ones(N_b,1);

%             A      = [sqrt(lambda_historyb)*eye(N_b,N_b); f_dyn; whatIsD];
            
            thisBehaviorSample = behavior{:,:}; 
            
            A = [zeros(N_b,N2) sqrt(lambda_historyb)*eye(N_b,N_b); sqrt(lambda_history)*eye(N2,N2) sqrt(lambda_history)*-1*f_dyn;whatIsD zeros(size(whatIsD,1),N_b);zeros(size(whatIsPsi,1),N2) sqrt(lambda_behavior)*whatIsPsi];
            Y = [sqrt(lambda_historyb)*bcoef_dcs(:,:,kk-1); zeros(N2,1); yNow; thisBehaviorSample(:,kk-1)];
            lamVec = [lambda_val*ones(N2,1);lambda_b*ones(N_b,1)];
            
    
            opts.pos        = false;
            opts.lambda     = lamVec(:);
            opts.check_grad = 0;
            res = fista_lasso(Y, A, [], opts);
    %         res = res + sign(res).*opts.lambda.*(abs(res) > 0.1*max(abs(res)));% Correct lasso bias
            % OPTION 1
    %         idxNNZ = (abs(res) > 0.5*max(abs(res)));
    %         res(~idxNNZ) = 0;
    %         res(idxNNZ)  = A(:,idxNNZ)\Y;
            % OPTION 2
    %         opts.lambda = lamVec(:).*(~idxNNZ);
            if debiasOpt
                opts.lambda = lamVec(:)./(1 + 100*abs(res));
                res         = fista_lasso(Y, A, [], opts);
            end
        end

    else % not Behavior

        
        if lambda_historyb == 0
            A = [zeros(N_b,N2) eye(N_b,N_b); sqrt(lambda_history)*eye(N2,N2) sqrt(lambda_history)*-1*f_dyn;whatIsD zeros(size(whatIsD,1),N_b)];
            Y = [bcoef_dcs(:,:,kk-1); zeros(N2,1); yNow];
            lamVec = [lambda_val*ones(N2,1);lambda_b*ones(N_b,1)];

            opts.pos        = false;
%             opts.lambda     = lambda_b*ones(N_b,1);
            opts.lambda     = lamVec(:);
            opts.check_grad = 0;

%             Y      = [bcoef_dcs(:,:,kk-1); zeros(N2,1); yNow];
% %             A      = [sqrt(lambda_historyb)*eye(N_b,N_b); f_dyn; whatIsD];
%             A = [zeros(N_b,N2) eye(N_b,N_b); eye(N2,N2) -1*f_dyn;whatIsD zeros(size(whatIsD,1),N_b)];
    
            res = fista_lasso(Y, A, [], opts); %FIXME: augmented Y needed, not just measurement
       
    %         A = linop_handles([M,M],Af,Ab,opt_set);
    %         res = solver_L1RLS(f_dyn, yNow,lambda_b,zeros(N_b,1),opts);
        else
            
%             A      = [sqrt(lambda_historyb)*eye(N_b,N_b); f_dyn];
%             Y      = [sqrt(lambda_historyb)*bcoef_dcs(:,:,kk-1); yNow];
%             lamVec = lambda_b*ones(N_b,1);

%             A      = [sqrt(lambda_historyb)*eye(N_b,N_b); f_dyn; whatIsD];
            A = [zeros(N_b,N2) sqrt(lambda_historyb)*eye(N_b,N_b); sqrt(lambda_history)*eye(N2,N2) sqrt(lambda_history)*-1*f_dyn;whatIsD zeros(size(whatIsD,1),N_b)];
            Y = [sqrt(lambda_historyb)*bcoef_dcs(:,:,kk-1); zeros(N2,1); yNow];
            lamVec = [lambda_val*ones(N2,1);lambda_b*ones(N_b,1)];
            
    
            opts.pos        = false;
            opts.lambda     = lamVec(:);
            opts.check_grad = 0;
            res = fista_lasso(Y, A, [], opts);
    %         res = res + sign(res).*opts.lambda.*(abs(res) > 0.1*max(abs(res)));% Correct lasso bias
            % OPTION 1
    %         idxNNZ = (abs(res) > 0.5*max(abs(res)));
    %         res(~idxNNZ) = 0;
    %         res(idxNNZ)  = A(:,idxNNZ)\Y;
            % OPTION 2
    %         opts.lambda = lamVec(:).*(~idxNNZ);
            if debiasOpt
                opts.lambda = lamVec(:)./(1 + 100*abs(res));
                res         = fista_lasso(Y, A, [], opts);
            end
        end

    end



    
    im_res              = DWT_invert(real(res(1:N2)));
    coef_dcs(:, :, kk)  = res(1:N2); 
    bcoef_dcs(:, :, kk) = res(N2+1:N2+N_b);


        
        % end % end fista  
  

    recon_dcs(:, :, kk) = im_res;
    if rMSE_calc_opt == 1
        rMSE_dcs(kk) = sum(sum((recon_dcs(:, :, kk) - TRUE_VID(:, :, kk)).^2))/sum(sum(TRUE_VID(:, :, kk).^2));
        PSNR_dcs(kk) = psnr(real(recon_dcs(:, :, kk)), TRUE_VID(:, :, kk), 2);
        TIME_ITER = toc;
        if verbose_flag == 1
            fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', kk, num_frames, TIME_ITER, PSNR_dcs(kk), rMSE_dcs(kk))
        end
    else
        TIME_ITER = toc;
        if verbose_flag == 1
            fprintf('Finished frame %d of %d in %f seconds.\n', kk, num_frames, TIME_ITER)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs

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

