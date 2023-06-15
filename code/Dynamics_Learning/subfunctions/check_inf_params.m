function inf_opts = check_inf_params(inf_opts)

% inf_opts = check_inf_params(inf_opts)
%
%
%
% 2018 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(inf_opts)                                                      % Make sure that bg_params is a struct
    clear inf_opts
    inf_opts = struct;
end
if (~isfield(inf_opts,'lambda_val'))||isempty(inf_opts.lambda_val)
    inf_opts.lambda_val = 0.8;
end
if (~isfield(inf_opts,'lambda_history'))||isempty(inf_opts.lambda_history)
    inf_opts.lambda_history = 0.9;
end
if (~isfield(inf_opts,'lambda_b'))||isempty(inf_opts.lambda_b)
    inf_opts.lambda_b = 0.02;
end
if (~isfield(inf_opts,'lambda_historyb'))||isempty(inf_opts.lambda_historyb)
    inf_opts.lambda_historyb = 0;
end
if (~isfield(inf_opts,'tol'))||isempty(inf_opts.tol)
    inf_opts.tol = 1e-3;
end
if (~isfield(inf_opts,'max_iters'))||isempty(inf_opts.max_iters)
    inf_opts.max_iters = 3e3;
end
if (~isfield(inf_opts,'plot_option'))||isempty(inf_opts.plot_option)
    inf_opts.plot_option = 100;
end
if (~isfield(inf_opts,'step_decay'))||isempty(inf_opts.step_decay)
    inf_opts.step_decay = 0.9999;
end
if (~isfield(inf_opts,'N_ex'))||isempty(inf_opts.N_ex)
    inf_opts.N_ex = 10;
end
if (~isfield(inf_opts,'step_f'))||isempty(inf_opts.step_f)
    inf_opts.step_f = 20;
end
if (~isfield(inf_opts,'step_d'))||isempty(inf_opts.step_d)
    inf_opts.step_d = 10;
end
if (~isfield(inf_opts,'T_s'))||isempty(inf_opts.T_s)
    inf_opts.T_s = 10;
end
if (~isfield(inf_opts,'ssim_flag'))||isempty(inf_opts.ssim_flag)
    inf_opts.ssim_flag = false;
end
if (~isfield(inf_opts,'dsamp'))||isempty(inf_opts.dsamp)
    inf_opts.dsamp = 4;
end
if (~isfield(inf_opts,'F_update'))||isempty(inf_opts.F_update)
    inf_opts.F_update = true;
end
if (~isfield(inf_opts,'lambda_f'))||isempty(inf_opts.lambda_f)
    inf_opts.lambda_f = 0;
end
if (~isfield(inf_opts,'lambda_f_decay'))||isempty(inf_opts.lambda_f_decay)
    inf_opts.lambda_f_decay = 1;
end
% if (~isfield(inf_opts,'solver_type'))||isempty(inf_opts.solver_type)
%     inf_opts.solver_type = 'tfocs';
% end
if (~isfield(inf_opts,'CVX_Precision'))||isempty(inf_opts.CVX_Precision)
    inf_opts.CVX_Precision = 'default'; % http://cvxr.com/cvx/doc/solver.html
end
if (~isfield(inf_opts,'deltaDynamics'))||isempty(inf_opts.deltaDynamics)
    inf_opts.deltaDynamics = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%