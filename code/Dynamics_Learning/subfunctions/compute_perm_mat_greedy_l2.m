function [ P, p_order ] = compute_perm_mat_greedy_l2( bhat, coeff_gt )
%COMPUTE_PERM_MAT_GREEDY [ P ] = compute_perm_mat_greedy( bhat, coeff_gt )
%   Greedy method to compute P that minimizes ||coeff_gt - P*bhat||_F^2,
%   such that |P(:,n)| = 1 and |P(n,:)| = 1


% b_indx corresponds to column of P, c_indx corresponds to row of P
eps = 1e-6;
[b_sz, T] = size(bhat);
c_sz = size(coeff_gt,1);
p_order = zeros(3,c_sz);
b_indx = 1:b_sz; c_indx = 1:c_sz;
P = zeros(c_sz,b_sz);
[A_mat,mat_err] = form_mats(bhat,coeff_gt);
mat_err0 = mat_err;

for ii = 1:b_sz
    [a, bx, cx] = min_err(A_mat,mat_err);
    P(c_indx(cx),b_indx(bx)) = a;
    p_order(:,ii) = [c_indx(cx); b_indx(bx); sign(a)];
    mat_err(cx,:) = nan; mat_err(:,bx) = nan;
    b_indx(bx) = nan; c_indx(cx) = nan;
end
    
    function [as, bx, cx] = min_err(A_mat,mat_err)
        [cx,bx] = find(mat_err == min(mat_err(:)));
        cx = cx(1); bx = bx(1);
        as = A_mat(cx,bx);
    end

    function [A_mat, mat_err] = form_mats(b_mat,c_mat)
        bsz = size(b_mat,1); csz = size(c_mat,1);
        mat_err = zeros(csz,bsz);
        A_mat = zeros(csz,bsz);
        for kk = 1:bsz
            for jj = 1:csz
                br = b_mat(kk,:); cr = c_mat(jj,:);
                as = (cr*br')/((br*br')+eps);   % adding eps here to avoid nan in mat_err
                A_mat(jj,kk) = as;
                mat_err(jj,kk) = (norm(cr-as*br,2)^2)/(norm(cr,2)^2);
            end
        end
    end

end

