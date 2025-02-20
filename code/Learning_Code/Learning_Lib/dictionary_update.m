function dict_new = dictionary_update(x_im, dict_old, coef_vals, step_s, varargin)

% function dict_new = dictionary_update(x_im, dictionary_old, coef_vals,
% step_s, opts)
% 
% Takes a gradient step with respect to the sparsity inducing energy
% function.
% 
% Inputs:
%   x_im        - Data samples over which to average the gradient step
%   dict_old    - The previous dictionary (used to infer the coefficients)
%   coef_vals   - The inferred coefficients for x_im using dict_old
%   step_s      - The step size to take in the gradient direction
%   opts        - Options for the particular problem (outlined in
%                 learn_dictionary.m)
%
% Outputs:
%   dict_new    - The new dictionary after the gradient step
% 
% Last Modified 6/4/2010 - Adam Charles

if nargin > 4
    opts = varargin{1};
else
    opts.in_iter = size(x_im, 2);
    opts.GD_iters = 1;
    opts.grad_type = 'norm';
    opts.nneg_dict = 0;
    opts.verysparsesuspected = 0; % very sparse suspected is useful when results are bad due to zeros but it's impossible to know without ground truth
    % opts.lambda2   = 0.1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take a gradient step

if strcmp(opts.grad_type, 'norm')
    for index2 = 1:opts.GD_iters
        % Take a step in the negative gradient of the basis:
        % Minimizing the energy:
        % E = ||x-Da||_2^2 + lambda*||a||_1^2
        % Update The basis matrix
        if opts.verysparsesuspected
            dict_new = dict_old + (step_s/(norm(coef_vals)*opts.in_iter))*...
                (x_im - dict_old*coef_vals)*coef_vals';
            warning('dictionary_update verysparsesuspected setting activated');
        else
            dict_new = dict_old + (step_s/opts.in_iter)*...
                (x_im - dict_old*coef_vals)*coef_vals';
        end
        % disp('In dictionary_update')
        % disp(size(x_im))
        % disp(size(coef_vals))
        % disp(size((step_s/opts.in_iter)*...
        %     (x_im - dict_old*coef_vals)*coef_vals'))
        % disp(size(dict_old))
        % disp(size(dict_new))

        % This part is basically the same, only for the
        % hyperspectral, care needs to be taken to saturate at 0,
        % so that no negative reflectances are learned. 
        if opts.nneg_dict == 1
            dict_new(dict_new < 0) = 0;
        end

        % % Re-normalize the basis
        dict_new = dict_new*diag(1./(sqrt(sum(dict_new.^2))));
        % if opts.verysparsesuspected ~= 1
        %     dict_new = dict_new*diag(1./(sqrt(sum(dict_new.^2))));
        % end
    end     
elseif strcmp(opts.grad_type, 'forb')
    for index2 = 1:opts.GD_iters
        % Take a step in the negative gradient of the basis:
        % This time the Forbenious norm is used to reduce unused
        % basis elements. The energy function being minimized is
        % then:
        % E = ||x-Da||_2^2 + lambda*||a||_1^2 + ||D||_F^2

        % Update The basis matrix
        % dict_new = dict_old + (step_s)*(...
        %     (x_im - dict_old*coef_vals)*coef_vals' -...
        %     opts.lambda2*2*dict_old)*diag(1./(1+sum(coef_vals ~= 0, 2)));

        % dict_new = dict_old + (step_s/(norm(coef_vals)*opts.in_iter))*(...
        %     (x_im - dict_old*coef_vals)*coef_vals' -...
        %     opts.lambda2*2*dict_old);

        % dict_new = dict_old + (step_s)*(...
        %     (x_im - dict_old*coef_vals)*coef_vals' -...
        %     opts.lambda2*2*dict_old)*diag(1./(norm(coef_vals)));

        dict_new = dict_old + (step_s)*(...
            (x_im - dict_old*coef_vals)*coef_vals' -...
            opts.lambda2*2*dict_old)*diag(1./(1+sum(coef_vals ~= 0, 2)));


        % For some data sets, the basis needs to be non-neg as well
        if opts.nneg_dict == 1
            dict_new(dict_new < 0) = 0;
        else
            % Do nothing
        end
    end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
