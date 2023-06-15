function [A,B] = rwl1df_bilinear_handle(X,D,F,opts)


param_vals.lambda_val = opts.lambda_val;
param_vals.gamma_val = opts.gamma_val; 
param_vals.alpha_a = opts.alpha_a;
param_vals.beta_a = opts.beta_a; 
param_vals.xi_a = opts.xi_a; 
param_vals.alpha_b = opts.alpha_b; 
param_vals.beta_b = opts.beta_b; 

if isfield(opts,'tol')
    param_vals.tol = opts.tol;
else
    param_vals.tol = 1e-3;
end

% f_dyn = @(z) F*z;
% 
% DYN_FUN{1} = f_dyn;

DWTfunc.apply = @(q) q;
DWTfunc.invert = @(q) q;

meas_func.Phi = @(z) D*z;
meas_func.Phit = @(z) (D')*z;
MEAS_FUN{1} = meas_func;

[A, ~] = RWL1_DF_bilinear(X, MEAS_FUN, F, DWTfunc, param_vals);

B = A{2};
A = A{1};

end
