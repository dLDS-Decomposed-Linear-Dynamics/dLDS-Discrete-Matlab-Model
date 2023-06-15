function varargout = rand_dyn_create(varargin)

% varargout = rand_dyn_create(varargin)
% 
% Creates a random set of dynamics in cell form
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 0
    f_size = varargin{1};
else
    error('You have to at least know what size vector you want!')
end

if nargin > 1
    f_num = varargin{2};
else
    f_num = 1;
end

if nargin > 2
    f_type = varargin{3};
else
    f_type = 'perm';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creat dynamics

F_cell = cell(1,f_num);                                                    % Initialize the cell array for the dynamics
for kk = 1:f_num
    if strcmp(f_type, 'perm') == 1
        F_cell{kk} = sparse(1:f_size,randperm(f_size), ...
                           rand_bern(f_size,0.5) + 0.05*randn(f_size,1));  % Create a random permutation matrix
        F_cell{kk} = full(F_cell{kk});                                     % Convert the matrix to a full matrix
    else
        error('unknown dynamics function type')
    end
end
F_cell = normalize_dynamics(F_cell);                                       % Normalize the dynamics
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout > 0
    varargout{1} = F_cell;
end

if nargout > 1
    for kk = 2:nargout
        varargout{kk} = [];
    end
else                                                                       % Nothing to do
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%