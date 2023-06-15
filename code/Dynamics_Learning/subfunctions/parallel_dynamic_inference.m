function A_cell = parallel_dynamic_inference(X, D, F, infer_hand, inf_opts, varargin)


for kk = 1:numel(X)
    X{kk} = reshape(X{kk},size(X{kk},1),1,size(X{kk},2)); 
end

if nargin < 6
    b = ones(1,numel(F));
else
    b = varargin{1};
end

F_mat = 0;
for kk = 1:numel(b)
    F_mat = F_mat + b(kk)*F{kk};
end

% Initialize outputs
A_cell = cell(size(X));

parfor kk = 1:numel(X)
    A_cell{kk} = feval(infer_hand, X{kk}, D, F_mat, inf_opts);
end

for kk = 1:numel(X)
    A_cell{kk} = reshape(A_cell{kk},size(A_cell{kk},1),size(A_cell{kk},3)); 
end

end
