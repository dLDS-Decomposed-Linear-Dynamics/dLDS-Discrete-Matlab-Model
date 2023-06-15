function [bcoef, varargout] = opt_regularized_AC(F, param_a, param_b, opts)

% This optimization attempts to optimize (minimize) the following cost
%
%              J(x) = ||A*x||_1 + sum_i log(|<u_i, b>| + nu)
%

% cvx_begin
% variable x(n);
% minimize(norm(A*x,1) - sum(log(abs(U_opt*x) + eta);
% subject to
%     U_opt*x > 0;
% cvx_end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Errors

if isfield(opts, 'eta')
    eta = opts.eta;
else
    error('Have to give an eta value!')
end

if isfield(opts, 'w1')
    w1 = opts.w1;
else
    w1 = 1;
end

if isfield(opts, 'w2')
    w2 = opts.w2;
else
    w2 = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alternative

% Create A
A = [w1*diag(param_a)*F; w2*diag(param_b)];
n = size(A,2);

% cvx_begin
%     variable x(n);
%     minimize(norm(A*x,1));
%     subject to
%         norm(x) <= 1;
% cvx_end

x0 = zeros(n,1);

norm_val = 2 - eta/mean(sqrt(sum(F.^2,2)));
n_samps = size(A,1)*2^(size(A,2));       % Choose this better later
x_samp = x0(:)*ones(1,n_samps) + randn(size(x0(:),1), n_samps);
x_samp = x_samp*sparse(1:n_samps, 1:n_samps, 1./sqrt(sum(x_samp.^2,1)))*norm_val;

eval_samp = sum(abs(A*x_samp),1) - sum(log(abs(F*x_samp) + eta));


opt_samp = x_samp(:, eval_samp == min(eval_samp));

if size(opt_samp, 2) > 1
%     opt_samp = opt_samp(:, ceil(rand(1)*size(opt_samp, 2)));
    opt_samp = opt_samp(:, 1);
end

% size(opt_samp)
% size(diag(sign(F*opt_samp)))

opt_cond = diag(sign(F*opt_samp))*F;

%cvx_begin
%    variable x(n);
%    minimize(norm(A*x,1) - sum(log(opt_cond*x + eta)) );
%    subject to
%        opt_cond*x >= 0;
%cvx_end
options = optimoptions('fmincon','Algorithm','interior-point');
x = fmincon(@(z) (sum(abs(A*z)) - sum(log(opt_cond*z+eta))), opt_samp, ...
    -opt_cond, zeros(size(opt_cond, 1), 1), [], [], [], [], [], options);

bcoef = x;

% if sum(isnan(bcoef))>0
%     sum(eval_samp == min(eval_samp))
% end

if nargout > 1
    varargout{1} = x_samp;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OLD
% This optimizatio takes 3 steps
% Step 1: Locate all the potential cones to optimize over
% Step 2: Assess which cone yields a good result
% Step 3: Perform a convex optimization of the cost function restricted to that cone


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Find all cones

% cone_set = find_nonint_cones(F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3
