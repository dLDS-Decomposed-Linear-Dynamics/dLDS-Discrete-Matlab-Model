function A = bpdndf_handle(X,D,F,opts)


param_vals.lambda_val = opts.lambda_val;   % 0.01;
param_vals.lambda_history = opts.lambda_history; % 0.3;
if isfield(opts,'tol')
    param_vals.tol = opts.tol;
else
    param_vals.tol = 1e-3;
end

f_dyn = @(z) F*z;

DYN_FUN{1} = f_dyn;

DWTfunc.apply = @(q) q;
DWTfunc.invert = @(q) q;

meas_func.Phi = @(z) D*z;
meas_func.Phit = @(z) (D')*z;
MEAS_FUN{1} = meas_func;

[A, ~] = BPDN_DF_largescale(X, MEAS_FUN, DYN_FUN, DWTfunc, param_vals);

end
