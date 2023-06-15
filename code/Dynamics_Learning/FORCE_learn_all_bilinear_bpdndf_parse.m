clc
clear

load('synth_rerun_FORCE_obs')

[data, F_out, coeff_gt] = rand_seq_create(sig_opts, noise_opts, F_true);

[ dhat, a, b ] = dL_reconstruct( data, D, F, infer_hand, inf_opts );

[P, p_order] = compute_perm_mat_greedy_l2(b(:,2:end), coeff_gt(:,2:end));
bt = P*b;
%%
figure(1)
clf
subplot(2,1,1)
plot(data')
subplot(2,1,2)
plot(dhat')

figure(2)
clf
plot(a')
title('obs coeff')

figure(3)
clf
subplot(3,1,1)
plot(coeff_gt')
title('dyn coeff gt')
subplot(3,1,2)
plot(b')
title('dyn coeff inf')
subplot(3,1,3)
plot((P*b)')
title('dyn coeff inf perm')

p_indx = 3;
figure(4)
clf
hold on
plot(coeff_gt(p_order(1,p_indx),2:end))
plot(bt(p_order(1,p_indx),2:end))
% plot(b(p_order(2,p_indx),2:end))
hold off
legend('coeff gt','b perm','b')
%%
