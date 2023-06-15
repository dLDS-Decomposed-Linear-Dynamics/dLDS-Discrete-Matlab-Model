function F = dynamics_update_for_DI(F, A_cell, step_F, B_cell,bias_vec,order)
% F_new = dynamics_update(F, A_cell, step_F, B_cell)
%
% Function to update the dyanmics operator dictionary F based on estiamted
% signal coefficents (A_cell) and dynamics coefficients (B_cell). The 
% update is via a gradient step (step size dictated by step_F) and includes
% a normalization sep to ensure that each operator is BIBO stable as a 
% discrete linear time-invariant system. 
%
% 2017 - Adam Charles

% Notes (Noga):
% 1)  B_cell  =  dynamics coefficients (varargin)
% 2) A_cell = signal coefficents 
% 3) A_now = initialization of signal coefficient array x_(t+1). It is
% actually X_n from y_n = \phi x_n +epsilon
% 4) E_all = matrix of error values (A_n - F_n * B_n)
% 5) A_del = delayed coefficients. It is x_(t)
% 6) F_delta = cell aray of operator gradients - gradient step for each dynamics operator 
% 7) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Input parsing 
%display(class(B_cell))
B_cell = B_cell{1};
A_cell = A_cell{1};
%display(class(A_cell))
%display(length(A_cell))
if nargin > 3  && ~isempty(B_cell)                                         % Check that the dynamics coefficients are provided
    % B_cell = B_cell;                                                     % If so directly read them in
else                                                                       % If not...
    B_cell = cell(1, numel(A_cell));                                       %   ... create a set of coefficients the correct size ...
    for kk = 1:numel(A_cell)                                               % 
        B_cell{kk} = ones(numel(F), size(A_cell{kk},2));                   %   ... and make all dynamics coefficients "1".
    end    
end
if nargin < 5 || isempty(bias_vec) || all(bias_vec == 0)
    bias_vec = zeros(size(A_cell,1),1) ;
end

x_vec = create_reco(F, B_cell, A_cell(:,1), 'ada', A_cell, 0,bias_vec,order,{},false);
x_vec = mean(cat(3,x_vec{:}),3);
disp(class(x_vec))
for f_spec_num = 1:length(F)
    grad_f_lose = zeros(size(F{1}));
    for time_point = 1:size(A_cell,2)
        grad_F_loss_spec = (A_cell(:,time_point+1) - x_vec(:,time_point+1)) * B_cell(f_spec_num,time_point) * x_vec(:,time_point)';
        grad_f_lose = grad_f_lose + grad_F_loss_spec;
    end
    grad_f_lose = grad_f_lose / size(A_cell,2);
    F{f_spec_num} = F{f_spec_num} - step_F * grad_f_lose;
    F_new = F{f_spec_num};
    max_svd = max(abs(svd(F_new)));                                        % Get the maximum singular value for the new dynamics function
    if max_svd == 0 %isnan(max_svd)||isinf(max_svd)                        % if the max singular value is 0...
        fprintf('WARNING: SVD had the value %f\n', max_svd);               %   ... and give a warning indicating a zero operator. 
    else
        F_new = F_new/max(abs(svd(F_new)));                    % Otherwise normalize by the operator norm 
    end                                                                    %
    F{f_spec_num} = F_new;
end                                                                        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Initializations
% 
% A_now = [];                                                                % Initialize the coefficient array
% A_del = [];                                                                % Initialize the delayedcoefficient array
% B_all = [];                                                                % Initialize the coefficient array
% % for kk = 1:numel(A_cell)                                                   % Iterate ove all example sequences (cells in the cell arrays)
% %     A_now = [A_now, A_cell{kk}(:, 2:end)];                                 % Add in the coefficients for each sequence
% %     A_del = [A_del, A_cell{kk}(:, 1:(end-1))];                             % Add in the delayed coefficients for each sequence
% %     B_all = [B_all, B_cell{kk}(:, 2:end)];                                 % Add in the dynamics coefficients for each sequence
% % end
% % [NM] change
% for kk = 1:numel(A_cell)                                                   % Iterate ove all example sequences (cells in the cell arrays)
%     A_now = [A_now, A_cell{kk}(:,2:size(B_cell{1},2))];                                 % Add in the coefficients for each sequence
%     A_del = [A_del, A_cell{kk}(:,1:size(B_cell{1},2)-1)];                             % Add in the delayed coefficients for each sequence
%     B_all = [B_all, B_cell{kk}];                                 % Add in the dynamics coefficients for each sequence
% end
% %disp(size(A_now))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Compute necessary values
% %disp('size A_cell');
% %disp(size(A_now));
% %disp(size(A_cell{kk}))
% 
% %disp(size(B_all))
% %disp(size(A_now))
% E_all = zeros(size(A_now));                                                % Initialize the matrix of error values
% for kk = 1:size(A_now,2)
%     %display('update_B_all_and_F');
%     %display(size(B_all(:, kk)));
%     %display(size(F));
%     F_tmp = dyn_cell2mat(F, B_all(:, kk));                                 % Make the matrix of all dynamics operators weighted by the dynamics coefficients
%     % [NM] - calculate the error in a DISCRETE dynamical system. The
%     % approximation is f*x_t, where the real values are x_(t+1)
%     
%     if ~include_bias
%         E_all(:, kk) = A_now(:, kk) - F_tmp*A_del(:, kk);                      % The error per example is the difference between the coefficients and coefficients predicted by the previous coefficients pushed forward via the dynamics
%     else
%         E_all(:, kk) = A_now(:, kk) - F_tmp*A_del(:, kk) - bias_vec;  
%     end
% end
% 
% F_delta = cell(size(F));                                                   % Initialize the cell aray of operator gradients
% for kk = 1:numel(F_delta)                                                  % Iterate over all elements of the dynamics dictionary
%     F_delta{kk} = 0;                                                       % Initialize the gradient step for each dynamics operator to 0.
% end
% 
% for kk = 1:size(E_all, 2)                                                  % Iterate over all the different data examples
%     tmp_mat = (E_all(:, kk)*(A_del(:, kk).'));                             % The derivative is the error multiplied by the coefficients (obtained by taking the gradient of least-squares terms)
%     for ll = 1:numel(F)                                                    % Iterate over all dynamics dictionary elements
%         % update delta F by adding it the coefficient of the dynamics multiplied
%         % by the gradient. It is in order to calculate the GD afterwards. 
%         %filling F_delta by the coeffieicient of each sub-dynamic
%         %multiplied by the  the error term
% 
%         % [NM} - here is the term of \nabla of the GD
%         F_delta{ll} = F_delta{ll} + B_all(ll, kk)*tmp_mat;                 % augment the 
%     end
% end
% % [NM] - below is the update by GD
% F_new = cell(size(F));                                                     % Initialize the new dynamics 
% for kk = 1:numel(F)                                                        % Iterate over all dynamics dictioanry elements
% 
%     F_new{kk} = F{kk} + (step_F/size(E_all,2))*F_delta{kk};                % Take a gradient step over each dynamics dictionary element (for GD)
%     F_new_kk_rep_na = F_new{kk};
%     % [NM] ADDITION
%     F_new_kk_rep_na(isnan(F_new_kk_rep_na)) = 0;
%     F_new{kk} = F_new_kk_rep_na;
%     
%     max_svd = max(abs(svd(F_new{kk})));                                    % Get the maximum singular value for the new dynamics function
%     if max_svd == 0 %isnan(max_svd)||isinf(max_svd)                        % if the max singular value is 0...
%         F_new{kk} = F{kk};                                                 %   ... just keep the dynamics function as all zeros ...
%         fprintf('WARNING: SVD had the value %f\n', max_svd);               %   ... and give a warning indicating a zero operator. 
%     else
%         F_new{kk} = F_new{kk}/max(abs(svd(F_new{kk})));                    % Otherwise normalize by the operator norm 
%     end                                                                    %
% end                                                                        %
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Set outputs
% 
% varargout{1} = F_new;
% varargout{2} = F_delta;
% 
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
