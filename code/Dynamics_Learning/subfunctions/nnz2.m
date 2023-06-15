function [ num ] = nnz2( v, tol )
%NNZ2 Summary of this function goes here
%   Detailed explanation goes here

num = numel(v) - sum(abs(v)<tol);

end

