function [x_out, time_tot] = modCS_small(varargin)

% [x_out, time_tot] = modCS_small(y, G, F, 0.5*tau, obs_var/dyn_var, 
%                                       maxiter); 
% 
% This function solves the BPDN optimization function on a time-varying
% signal under the Gaussian innovations and sparse state assumption. Under
% these assumptions, BPDN is solved for the state at every time-step
% with the previous state used as additional measurements.
% 
%   The inputs are:
% 
% y:       MxT matrix of measurement vectors at each time step
% G:       MxNxT array of measurement matrices at each iteration
% F:       MxMxT array of dynamics matrices for each iteration
% tau:     1x2 vector containing the sparsity tradeoff parameter for BPDN 
%          wrt the innovations and the ratio of the innovations variance to
%          the measurement noise variance
% maxiter: Maximum number of iterations for BPDN to run
% 
%    The outputs are:
% x_out:      NxT matrix with the estimates of the signal x
% time_tot:   total runtime
% 
%
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated August 21, 2012. 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing
if nargin > 1
    y = varargin{1};
    G = varargin{2};
else
    error('Bad number of inputs!')
end
if nargin > 2
    F = varargin{3};
else
    F = eye(size(G, 2));
end

if nargin > 3
    tau = varargin{4};
else
    tau = 0.01;
end

if nargin > 4
    std_quot = varargin{5};
else
    std_quot = 1;
end

if nargin > 5
    maxiter = varargin{6};
else
    maxiter = 1000;
end

if nargin > 6
    error('Bad number of inputs!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations
x_dim = size(F, 1);
y_dim = size(y, 1);
num_iters = size(y, 2) - 1;

opts.tol = 0.001;
opts.printEvery = 0;

x_out = zeros(x_dim, num_iters+1);
time_bp = zeros(num_iters + 1);

% Initial Point
[x_out(:, 1), ~, ~, time_bp(1)] = BPDN_homotopy_function(G(:, :, 1), y(:, 1), std_quot, maxiter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Algorithm
% Iterate through
for ii = 1:num_iters
    % Optimize the BPDN-DF objective
    Af   = @(x) G(:, :, ii+1)*x;
    Ab   = @(x) (G(:, :, ii+1).')*x;
    A    = linop_handles([y_dim, x_dim], Af, Ab, 'R2R');
    
%     x_pred  = DWT_apply(f_dyn(vid_recon_dcs(:, :, kk-1)));
    weights = abs(F(:, :, ii+1)*x_out(:, ii))>=0.01;
    Wapply  = @(z) weights(:).*z(:);
    Wtrans  = @(z) weights(:).*z(:);
    W       = linop_handles([x_dim,x_dim], Wtrans, Wapply );
    
    % fprintf('Sparsity = %f\n', sum(weights)/numel(weights));
    
    mu = 0.1;
    x0 = [];
    z0 = [];
    opts.continuation = true;
    x_out(:, ii+1) = solver_sBPDN_W(A, W,...
      y(:, ii+1), tau, mu, x0, z0, opts );
end

time_tot = sum(time_bp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
