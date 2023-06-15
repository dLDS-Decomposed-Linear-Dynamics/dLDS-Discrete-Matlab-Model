function [x_supp, G_supp, F_supp] = restrict_model(x, G, F)

% [x_supp, G_supp, F_supp] = restrict_model(x, G, F)
% 
% Restrics the sparse dymic model to a model where the support of each
% signal is known. 
% 
% 
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated October 12, 2012. 
% 

M = size(G, 1);
P = size(x, 2);
S = max(sum(x ~= 0, 1));

x_supp = zeros(S, P);
G_supp = zeros(M, S, P);
F_supp = zeros(S, S, P);

for kk = 1:P-1
    supp_loc = x(:, kk) ~= 0;
    Sn = sum(supp_loc);
    x_supp(:, kk) = [x(x(:, kk) ~= 0, kk); zeros(S-Sn, 1)];
    G_supp(:, :, kk) = [G(:, x(:, kk) ~= 0, kk), zeros(M, S-Sn)];
    if kk > 1
        F_supp(:, :, kk) = [F(x(:, kk) ~= 0, x(:, kk-1) ~= 0, kk), zeros(Sn, S - Sn_old); zeros(S - Sn, S)];
    end
    Sn_old = Sn;
end

end