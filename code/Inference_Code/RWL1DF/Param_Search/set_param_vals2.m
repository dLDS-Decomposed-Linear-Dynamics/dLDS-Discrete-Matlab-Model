function param_vals = set_param_vals(TOL, mult_range, den_range)


if nargin > 3

    if nargin > 4
        param_vals.beta = varargin{1};
    else
        param_vals.beta = 1;
    end
    param_vals.tol = TOL;


    param_vals.rwl1_mult =  mult_range;
    param_vals.rwl1_reg = den_range;
    param_vals.lambda_val = lambda_range;

elseif nargin <= 3

     param_vals.lambda_val = mult_range;
     param_vals.EPS = den_range;
     param_vals.tol = TOL;

else

% Do nothing

end


end
