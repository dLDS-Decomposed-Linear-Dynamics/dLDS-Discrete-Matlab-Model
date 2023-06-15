function min_params = find_opt_params(RES_MAT, varargin)


dim_num = numel(size(RES_MAT));

if nargin < dim_num
	error('Bad numbers of inputs!')
end

param_ranges = cell(dim_num-1,1);
for kk = 2:dim_num
	param_ranges{kk-1} = varargin{kk-1};
end

mean_mat = mean(RES_MAT, dim_num);
min_RES = min(mean_mat(:));

res_szs = size(mean_mat);
min_ind = find(min_RES == mean_mat);

ind_str = '[I1';
for kk = 1:(dim_num-1)
	ind_str = [ind_str,sprintf(',I%d',kk+1)];
end
ind_str = [ind_str,']'];
eval_str = [ind_str,' = ind2sub(res_szs,min_ind);'];
eval(eval_str);

eval(['param_locs = ',ind_str,';']);

for kk = 1:(dim_num-1)
	min_params(kk) = param_ranges{kk}(param_locs(kk));
end

end
