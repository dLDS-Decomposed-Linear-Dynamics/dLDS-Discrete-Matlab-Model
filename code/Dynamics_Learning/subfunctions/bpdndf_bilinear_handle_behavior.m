function [A,B] = bpdndf_bilinear_handle_behavior(X,D,F,opts, varargin)


% param_vals.lambda_val      = opts.lambda_val;                              % extract lambda value;
% param_vals.lambda_b        = opts.lambda_b;                                % extract dynsmics-lambda value;
% param_vals.lambda_history  = opts.lambda_history;                          % extract history lambda value;
% param_vals.lambda_historyb = opts.lambda_historyb;                         % extract dynamics histroy lambda value;

param_vals = opts;

if ~isfield(param_vals,'tol')
    param_vals.tol = 1e-3;
end

if nargin > 4
    Psi = varargin{1};
    behavior = varargin{2};
end

% f_dyn = @(z) F*z;
% 
% DYN_FUN{1} = f_dyn;

DWTfunc.apply = @(q) q;
DWTfunc.invert = @(q) q;

meas_func.Phi  = @(z) D*z;
meas_func.Phit = @(z) (D')*z;
meas_func.Psi  = @(z) Psi*z;
meas_func.Psit = @(z) (Psi')*z;

MEAS_FUN{1} = meas_func;

if strcmp(opts.special,'sysop')
    [A, ~] = BPDN_DF_largescale(X, MEAS_FUN, @(z) F{1}*z, DWTfunc, param_vals);
    B = ones(1,1,size(A,3));
elseif strcmp(opts.special,'switched')
    [A, ~] = BPDN_DF_switched(X, MEAS_FUN, F, DWTfunc, param_vals);
    B = A{2};
    A = A{1};
elseif strcmp(opts.special,'nofilt')
    [A, ~] = BPDN_DF_bilinearNoFilt(X, MEAS_FUN, F, DWTfunc, param_vals);
    B = A{2};
    A = A{1};
elseif strcmp(opts.special,'noobs')
    [A, ~] = BPDN_DF_bilinearNoObs(X, MEAS_FUN, F, DWTfunc, param_vals);
    B = A{2};
    A = A{1};
elseif strcmp(opts.special,'nodyn')
    [A, ~] = BPDN_DF_bilinearNoDyn(X, MEAS_FUN, F, DWTfunc, param_vals);
    B = A{2};
    A = A{1};
elseif opts.behaviordLDS == 1
    if strcmp(opts.special,'noobs')
        [A, ~] = BPDN_DF_bilinearNoObs(X, MEAS_FUN, F, DWTfunc, param_vals,Psi,behavior);
        B = A{2};
        A = A{1};
    else
        [A, ~] = BPDN_DF_bilinear_behavior(X, MEAS_FUN, F, DWTfunc, param_vals,Psi,behavior);
        B = A{2};
        A = A{1};
    end
else
    [A, ~] = BPDN_DF_bilinear(X, MEAS_FUN, F, DWTfunc, param_vals);
    B = A{2};
    A = A{1};
end



end
