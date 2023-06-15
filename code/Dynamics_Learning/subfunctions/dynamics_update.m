function varargout = dynamics_update(F, A_cell, step_F, varargin)

% F_new = dynamics_update(F, A_cell, step_F, B_cell)
%
% Function to update the dyanmics operator dictionary F based on estiamted
% signal coefficents (A_cell) and dynamics coefficients (B_cell). The 
% update is via a gradient step (step size dictated by step_F) and includes
% a normalization sep to ensure that each operator is BIBO stable as a 
% discrete linear time-invariant system. 
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing 

if nargin > 3                                                              % Check that the dynamics coefficients are provided
    B_cell = varargin{1};                                                  % If so directly read them in
else                                                                       % If not...
    B_cell = cell(1, numel(A_cell));                                       %   ... create a set of coefficients the correct size ...
    for kk = 1:numel(A_cell)                                               % 
        B_cell{kk} = ones(numel(F), size(A_cell{kk},2));                   %   ... and make all dynamics coefficients "1".
    end    
end

if nargin > 4;    lambda_f = varargin{2};
else;             lambda_f = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

A_now = [];                                                                % Initialize the coefficient array
A_del = [];                                                                % Initialize the delayedcoefficient array
B_all = [];                                                                % Initialize the coefficient array
for kk = 1:numel(A_cell)                                                   % Iterate ove all example sequences (cells in the cell arrays)
    A_now = cat(2, A_now, A_cell{kk}(:, 2:end));                                 % Add in the coefficients for each sequence
    A_del = cat(2, A_del, A_cell{kk}(:, 1:(end-1)));                             % Add in the delayed coefficients for each sequence
    B_all = cat(2, B_all, B_cell{kk}(:, 2:end));                                 % Add in the dynamics coefficients for each sequence
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute necessary values

nGrad = 1;

for gg = 1:nGrad
    F_delta = cell(size(F));                                                   % Initialize the cell aray of operator gradients
    for kk = 1:numel(F_delta)                                                  % Iterate over all elements of the dynamics dictionary
        F_delta{kk} = 0;                                                       % Initialize the gradient step for each dynamics operator to 0.
    end

    % The error terms are:
    %    E       = J(A) + sum_t||At - (sum Fk*ck)A{t-1}||_2^2
    %    dE/dFk  = d/dFk sum_t ||(At - sum Fk*(ck*A{t-1})||_2^2
    %            = sum_t d/dFk -2*ck*(At - sum Fk*(ck*A{t-1}))*A{t-1}.'

    E_all = zeros(size(A_now));                                                % Initialize the matrix of error values
    for kk = 1:size(A_now,2)
        F_tmp = dyn_cell2mat(F, B_all(:, kk));                                 % Make the matrix of all dynamics operators weighted by the dynamics coefficients
        E_all(:, kk) = A_now(:, kk) - F_tmp*A_del(:, kk);                      % The error per example is the difference between the coefficients and coefficients predicted by the previous coefficients pushed forward via the dynamics
    end

    for kk = 1:size(E_all, 2)                                                  % Iterate over all the different data examples
        tmp_mat = (E_all(:, kk)*(A_del(:, kk).'));                             % The derivative is the error multiplied by the coefficients (obtained by taking the gradient of least-squares terms)
        for ll = 1:numel(F)                                                    % Iterate over all dynamics dictionary elements
            F_delta{ll} = F_delta{ll} + B_all(ll, kk)*tmp_mat;                 % augment the 
        end
    end

    F_new = cell(size(F));                                                     % Initialize the new dynamics 
    for kk = 1:numel(F)                                                        % Iterate over all dynamics dictioanry elements
        F_new{kk} = F{kk} + (step_F/size(E_all,2))*F_delta{kk};                % Take a gradient step over each dynamics dictionary element
        F_new{kk} = prjectDown(F_new{kk}, lambda_f, F, kk, 'corr');                                                                    %
    end                                                                        %
    F = F_new;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs

varargout{1} = F_new;
varargout{2} = F_delta;

end

%%

function F_new = softThreshold(F_new, lambda_f)

IX         = (abs(F_new) < lambda_f);
F_new(IX)  = 0;
F_new(~IX) = F_new(~IX) - sign(F_new(~IX))*lambda_f;

end

function F_new = hardThreshold(F_new, lambda_f)

F_new(abs(F_new) < lambda_f) = 0;

end


function F_new = prjectDown(F_new, lambda_f, F, kk, opt)

switch lower(opt)
    case 'forb'
        F_new = F_new - lambda_f*F{kk};
    case 'l1'
        F_new = softThreshold(F_new, lambda_f);
    case 'l0'
        F_new = hardThreshold(F_new, lambda_f);
    case 'nuc'    
    case 'corr'
%         F_new = F_new - 0.0001*F{kk};
        for ll = 1:numel(F)
            if ll~=kk
                F_new = F_new - lambda_f*trace((F{kk}.')*F{ll})*F{ll};
            end
        end
    otherwise
        % do nothing
end

max_svd = max(abs(svd(F_new)));                                            % Get the maximum singular value for the new dynamics function
if max_svd == 0 %isnan(max_svd)||isinf(max_svd)                            % if the max singular value is 0...
    F_new = F{kk};                                                         %   ... just keep the dynamics function as all zeros ...
    fprintf('WARNING: SVD had the value %f\n', max_svd);                   %   ... and give a warning indicating a zero operator. 
else                                                                       % Otherwise normalize by the operator norm
%     if lambda_f > 0
%         F_new = F_new./max(1, max_svd);                                    % if lambda_f exists, allow dynamics to go to zero
%     else
        F_new = F_new./max_svd;                                            % Otherwise just normalize by the operator norm 
%     end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
