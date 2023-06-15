function F_err = correlate_learned_dynamics(F_learned, F_true,varargin)

% F_err = correlate_learned_dynamics(F_learned, F_true)
% 
% Correlate the learned dynamcis with a set of known dynamics. 
% 
% 
% 
% 
% August 13, 2014
% Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing
if nargin < 2
    error('Too few inputs!');
end

if nargin < 3
    D = eye(size(F_true{1}, 1));
end
if nargin > 2
    D = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Calculations

corr_tmp = zeros(numel(F_learned), numel(F_true));
F_errs_tmp = zeros(numel(F_learned), numel(F_true));

for kk = 1:numel(F_learned)
    for ll = 1:numel(F_true)
        % Calculate inner product
        corr_tmp(kk,ll) = sum(sum((D*F_learned{kk}*D').*F_true{ll}));
        % Calculate error
        F_errs_tmp(kk,ll) = sum(sum((D*F_learned{kk}*D' - F_true{ll}).^2))./sum(sum(F_true{ll}.^2));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlate Dynamics

F_err = zeros(numel(F_learned), 1);
idx = 1:numel(F_learned);
idy = 1:numel(F_true);

% Find correlations
while numel(idx) > 0
    [idx_tmp, idy_tmp] = find(corr_tmp == max(corr_tmp(:)));
    idx_tmp = idx_tmp(1);
    idy_tmp = idy_tmp(1);
    F_err(idx(idx_tmp)) = F_errs_tmp(idx(idx_tmp),idy(idy_tmp));
    idx = idx([1:idx_tmp-1,idx_tmp+1:end]);
    idy = idy([1:idy_tmp-1,idy_tmp+1:end]);
    corr_tmp = corr_tmp([1:idx_tmp-1,idx_tmp+1:end], [1:idy_tmp-1,idy_tmp+1:end]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
