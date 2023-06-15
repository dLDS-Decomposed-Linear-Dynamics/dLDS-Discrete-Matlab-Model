function [ P, p_order ] = compute_perm_mat_greedy_corr( bhat0, coeff_gt0 )
%COMPUTE_PERM_MAT_GREEDY [ P ] = compute_perm_mat_greedy( bhat, coeff_gt )
%   Greedy method to compute P that minimizes ||coeff_gt - P*bhat||_F^2,
%   such that |P(:,n)| = 1 and |P(n,:)| = 1


bhat = bhat0; coeff_gt = coeff_gt0;

% b_indx corresponds to column of P, c_indx corresponds to row of P

[b_sz, T] = size(bhat);
c_sz = size(coeff_gt,1);
p_order = zeros(2,c_sz);
P = zeros(c_sz,b_sz);
for ii = 1:b_sz
    corr_tmp = bhat*coeff_gt';
    [bx, cx] = find( abs(corr_tmp) == max(abs(corr_tmp(:))));
    b = bhat(bx,:); c = coeff_gt(cx,:);
    a = (c*b')/(b*b');
    P(cx,bx) = a;
    bhat(bx,:) = nan(1,T); coeff_gt(cx,:) = nan(1,T);
    p_order(:,ii) = [cx; bx];
end

end

