
lam_sweep = 1.5:0.1:2;
t_sweep = 0.1:0.2:1;
% Pi0_sweep = 1:2:10;
% d_sweep = 0.1:0.2:1;


MSE_TMP = zeros(numel(lam_sweep), numel(t_sweep));
MAX_TMP = zeros(numel(lam_sweep), numel(t_sweep));
MIN_TMP = zeros(numel(lam_sweep), numel(t_sweep));

for kk = 1:numel(lam_sweep)
    for ll = 1:numel(t_sweep)
% for kk = 1:numel(Pi0_sweep)%numel(lam_sweep)
%     for ll = 1:numel(d_sweep)%numel(t_sweep)
%         Pi0 = Pi0_sweep(kk);
%         lambdap = 1.7; %lam_sweep(kk);
        Pi0 = 1;
        lambdap = lam_sweep(kk);
%         [x_lscs, ~] = kfcsls_full(y, F, Pi0, G, 0.5*eye(x_dim), (t_sweep(ll))*eye(y_dim),lambdap);
%         [x_lscs, ~] = kfcsls_full(y, F, Pi0, G, 0.5*eye(x_dim), (0.3)*eye(y_dim),lambdap);
        %RWL1_DFmult(y, G, F, 100.*(dyn_var), obs_var, 20, tau, 0.01);
        x_lscs_all(:, :, 1) = x_lscs;
        MSE_TMP(kk, ll) =  mean(mean(sum((x_store - x_lscs_all).^2, 1)./x_norms, 3).');
        MIN_TMP(kk, ll) =  min(mean(sum((x_store - x_lscs_all).^2, 1)./x_norms, 3).');
        MAX_TMP(kk, ll) =  max(mean(sum((x_store - x_lscs_all).^2, 1)./x_norms, 3).');
        fprintf('Ran lscs %d %d : min = %f, mean = %f, max = %f...\n', kk, ll, MIN_TMP(kk, ll), MSE_TMP(kk, ll), MAX_TMP(kk, ll))
    end
end