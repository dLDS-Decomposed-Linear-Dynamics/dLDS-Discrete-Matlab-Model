function [ P ] = gen_rand_perm_mat( rn, cn )
%GEN_RAND_PERM_MAT [ P ] = gen_rand_perm_mat( rn, cn )
%   Generate random permutation and rescaling matrix

sz = min(rn,cn);
P = zeros(rn,cn);

r_indx = 1:rn; c_indx = 1:cn;

for ii = 1:sz
    rr = randperm(numel(r_indx),1); cr = randperm(numel(c_indx),1);
    P(r_indx(rr),c_indx(cr)) = randn;
    r_indx(rr) = []; c_indx(cr) = [];
end

