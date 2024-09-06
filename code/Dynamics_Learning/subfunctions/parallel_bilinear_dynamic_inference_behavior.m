function [A_cell,B_cell] = parallel_bilinear_dynamic_inference_behavior(X, D, F, Psi, behavior, infer_hand, inf_opts, varargin)

% [A_cell,B_cell] = parallel_bilinear_dynamic_inference(X, D, F, infer_hand, inf_opts, varargin)
%
% Function to infer the representational and dynamics coefficients under
% the BPDN-DF model. 
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up initializations 

for kk = 1:numel(X)
    X{kk} = reshape(X{kk},[size(X{kk},1),1,size(X{kk},2)]); 
    behavior{kk} = reshape(behavior{kk},[size(behavior{kk},1),1,size(behavior{kk},2)]); 

end

A_cell = cell(size(X));                                                    % Initialize representational coefficients
B_cell = cell(size(X));                                                    % Initialize dynamics coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run inference

parfor kk = 1:numel(X)
%for kk = 1:numel(X)
    [A_cell{kk},B_cell{kk}] = feval(infer_hand, X{kk}, D, F, inf_opts, Psi, behavior);    % Perform inference over each sequence independently
end

for kk = 1:numel(X)
    A_cell{kk} = reshape(A_cell{kk},size(A_cell{kk},1),size(A_cell{kk},3));% Reshape representational coefficients
    B_cell{kk} = reshape(B_cell{kk},size(B_cell{kk},1),size(B_cell{kk},3));% Reshape dynamics coefficients
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%