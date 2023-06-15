function [ sup, sup_mat ] = support_vectorize( F, th )
%SUPPORT_VECTORIZE [ sup ] = support_vectorize( F, th )
%   Find support of F. Zero out any values below th

if nargin < 2
    th = 0;
end

nF = numel(F);
nF1 = numel(F{1});

sup_mat = zeros(nF1,nF);

for ii = 1:nF
    Ft = F{ii};
    Ft(abs(Ft) < th) = 0;
    s = Ft ~= 0;
    sup_mat(:,ii) = s(:);
end
sup = sup_mat(:);

