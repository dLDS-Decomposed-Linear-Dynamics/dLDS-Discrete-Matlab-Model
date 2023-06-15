function F_mat = dyn_cell2mat(F,b)

% F_mat = dyn_cell2mat(F,b)
%
% Function to create the overall dynamics matrix 
%
% 2019 - Adam Charles

F_mat = 0;                                                                 % Initialize the full function matrix to zero
for kk = 1:numel(F)                                                        % Loop over all dynamics functions (elements of the cell array F)
    F_mat = F_mat + F{kk}*b(kk);                                           % Add the contribution of the k^th dynamics
end

end