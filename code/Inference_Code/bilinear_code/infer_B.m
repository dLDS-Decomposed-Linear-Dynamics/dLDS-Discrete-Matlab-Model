function [tot_coeffs, tot_reco_dyn,fx_cell] = infer_B(base_F, type_solve, clean_dyn, inf_opts)
% type solve can be:
% inv - find c's by inverse
% lasso - find c's by lasso
% tfocs - using tfocs


if nargin < 2 || isempty(type_solve), type_solve = 'inv'; end
if nargin < 4 || isempty(inf_opts)
    inf_opts = struct ;
    inf_opts.reg_c = 0.02;
end
%if nargin < 5 || isempty(bias_in_reco), bias_in_reco = false; end
if length(size(clean_dyn)) == 3
    clean_dyn = permute(clean_dyn,[1,3,2]);
end
time_points = size(clean_dyn ,2);
tot_coeffs = zeros(length(base_F), time_points-1);
tot_reco_dyn = zeros(size(clean_dyn));
tot_reco_dyn(:,1) = clean_dyn(:,1);
fx_cell ={};
%reco_error = zeros(time_points-1);
for time_pont = 1:time_points-1
    cur_dyn = clean_dyn(:,time_pont);
    next_dyn = clean_dyn(:,time_pont+1);

    if isfield(inf_opts,'include_bias') && inf_opts.include_bias
       next_dyn = next_dyn - inf_opts.bias_vec; 
    end
    fx = [];
    for f_num = 1:length(base_F)
        f_cur = base_F{f_num};
        %cur_x = return_struct.x{1};
        fx(:,end+1) = f_cur*cur_dyn;        
    end

    if strcmp(type_solve,'inv')
        coeffs = pinv(fx)*next_dyn;
    elseif strcmpi(type_solve,'lasso')
        coeffs = lasso(fx,next_dyn,'Lambda',inf_opts.reg_c);

    elseif strcmpi(type_solve,'tfocs')
        
        [ coeffs, ~, ~ ] = solver_L1RLS( fx, next_dyn,inf_opts.reg_c );
    else
        error('Unknown solver type')
    end
    reco_dyn = fx*coeffs;
    %disp(fx*coeffs)
    if isfield(inf_opts,'include_bias') && inf_opts.include_bias && inf_opts.include_reco_bias 
        tot_reco_dyn(:,time_pont+1) = reco_dyn+inf_opts.bias_vec;
    else
        tot_reco_dyn(:,time_pont+1) = reco_dyn;
    end
    tot_coeffs(:,time_pont) = coeffs;
    fx_cell{time_pont} = fx;
    %reco_dyn =
    %reco_error(time_pont) = 
end
end



