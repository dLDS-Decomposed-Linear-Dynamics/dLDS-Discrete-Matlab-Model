function [ img_cell ] = pairwise_dict_proj_wrapper( F1, F2, D, frm, proj_opts )
%PAIRWISE_DICT_PROJ_WRAPPER [ img_cell ] = pairwise_dict_proj_wrapper( F1, F2, D, frm, T, theta, scalar, scalar_method )
%   Detailed explanation goes here

T = proj_opts.dyn_steps;
theta = proj_opts.dict_coeff;
scalar = proj_opts.coeff_scalar;
scalar_method = proj_opts.scalar_method;
reshape_dim = proj_opts.reshape_dim;

img_cell = cell(size(theta,2),T);

for ss = 1:size(theta,2)
    th = theta(:,ss);
    proj = pairwise_dict_proj( F1, F2, frm, T, th, scalar, scalar_method );
    img = D*proj;
    for tt = 1:T
        img_cell{ss,tt} = reshape(img(:,tt),reshape_dim);
    end
end

