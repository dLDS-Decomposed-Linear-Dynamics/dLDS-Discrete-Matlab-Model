function [F, b] = mk_bilinear_F(F_cell, f_s, varargin)

% [F, b] = mk_bilinear_F(F_cell, f_s, varargin)
% 
% Create a random set of sparse weights and then generate a linear sum of
% the dynamics function based on those weights.
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 2
    unit_b = varargin{1};
else
    unit_b = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the bilinear weights

if numel(F_cell) == 1
    b = 1;                                                                 % If only one fynamics function to choose from, just use that
else
    b = zeros(numel(F_cell), 1);                                           % Otherwise initialize to zeros
    b(randsample(numel(b),f_s)) = randn(f_s, 1);                           % Set a random subset of the weights to non-zero values
    if unit_b
        b(b ~=0) = 1;                                                      % If needed, set the weight of each 'b' chosen to one 
    else
        b = b/sqrt(sum(b.^2));                                             % Otherwise just normalize the weights
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the full dynamics matrix

F = 0;                                                                     % Initialize the dynamics to zero
for kk = 1:numel(F_cell)
    if b(kk)~= 0
        F = F + F_cell{kk}*b(kk);                                          % Iteratively add in each dynamcs matrix's contribution
    end
end	

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%