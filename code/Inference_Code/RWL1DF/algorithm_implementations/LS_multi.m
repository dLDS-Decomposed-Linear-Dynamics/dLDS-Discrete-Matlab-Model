function [x_out, time_tot] = LS_multi(varargin)

% [x_out, time_tot] = BPDN_multi(y, G, tau, 1000);
% 
% This function solves the BPDN optimization function on a time-varying
% signal by solving for each time-step independently.
% 
%   The inputs are:
% 
% y:       MxT matrix of measurement vectors at each time step
% G:       MxNxT array of measurement matrices at each iteration
% tau:     1x2 vector containing the sparsity tradeoff parameter for BPDN 
%          wrt the innovations and the ratio of the innovations variance to
%          the measurement noise variance
% maxiter: The maximum number of iterations for the BPDN solver
% 
%    The outputs are:
% 
% x_out:      NxT matrix with the estimates of the signal x
% time_tot:   total runtime
% 
%
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated August 14, 2012. 
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
    supp_res = varargin{3};
else
    supp_res = 1;
end

if nargin > 3
    error('Bad number of inputs!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations
x_dim = size(G, 2);
num_iters = size(y, 2) - 1;

x_out = zeros(x_dim, num_iters+1);
time_bp = zeros(num_iters + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve all LS optimizations
% Can do this in parallel with parfor
parfor ii = 1:num_iters+1
    % Calculate pseudo inverse
    if supp_res == 1
        tic
        G_norm = sum(abs(G(:, :, ii)), 1);
        G_num = sum(G_norm~=0);
        G_tmp = G(:, 1:G_num, ii);
        x_tmp = G_tmp\y(:, ii);
        x_out(:, ii) = [x_tmp; zeros(x_dim-G_num, 1)];
        time_bp(ii) = toc;
    else
        tic
        x_out(:, ii) = G(:, :, ii)\y(:, ii);
        time_bp(ii) = toc;
    end
end

time_tot = sum(time_bp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
