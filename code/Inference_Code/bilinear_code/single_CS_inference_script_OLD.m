%% Mini Inference script
clear
% load /home/adam/GITrepos/dynamics_learning/results/OLD/20150216_BBCvideo_12x12_4x_nF20.mat
load /home/adam/GITrepos/dynamics_learning/results/OLD/20150210_BBCvideo_8x8_4x_nF20.mat

x_sel = 10;

n_meas = ceil(0.2*size(D,1));
% n_meas = 12

X_ex = cell(1, x_sel);
% Get a random video
for ll = 1:x_sel
    while sum(sum(abs(X_ex{ll}))==0)>0
        X_ex{ll} = rand_bbc_video([sqrt(sig_opts.M), sqrt(sig_opts.M), sig_opts.T_s], 2);
    end
end

G = cell(1,numel(X_ex)); 
Z = cell(1,numel(X_ex));

for ll = 1:numel(X_ex)
%     Z = zeros(n_meas,1,sig_opts.T_s);
    for kk = 1:sig_opts.T_s
        % Make a random 
        % G = eye(size(D, 1));
        Gtmp = randn(n_meas, size(D, 1));
        Gtmp = Gtmp*diag(1./sqrt(sum(Gtmp.^2,1)));    
        G{ll}(:,:,kk) = Gtmp;
        Z{ll}(:,:,kk) = G{ll}(:,:,kk)*X_ex{ll}(:,kk) + 0.0005*randn(n_meas,1);
    end
end 

opts = inf_opts;
param_vals.lambda_val      = opts.lambda_val;   % 8*;
param_vals.lambda_b        = 0.01*opts.lambda_b;   % 4*;
param_vals.lambda_history  = opts.lambda_history; % 2*;
param_vals.lambda_historyb = opts.lambda_historyb; %
if isfield(opts,'tol')
    param_vals.tol = opts.tol;
else
    param_vals.tol = 1e-3;
end

MEAS_FUN = cell(1,sig_opts.T_s);
DWTfunc.apply = @(q) q;
DWTfunc.invert = @(q) q;

F_1 = cell(1);
F_1{1} = eye(size(D,2));
f_dyn = @(q) q;

A = cell(1, numel(X_ex));
A2 = cell(1, numel(X_ex));
err_d = zeros(numel(X_ex), sig_opts.T_s);
err_nd = zeros(numel(X_ex), sig_opts.T_s);

start_ss = 1;

%%

for ll = 1:numel(X_ex)
    
    for kk = 1:sig_opts.T_s
        meas_func.Phi = @(z) G{ll}(:,:,kk)*D*z;
        meas_func.Phit = @(z) ((G{ll}(:,:,kk)*D)')*z;
        MEAS_FUN{kk} = meas_func;
    end
    
    % Modify Parameters?
    param_vals.lambda_val      = 0.01*opts.lambda_val; %max(max((G{ll}(:,:,1)')*Z{ll}(:,:,1).^2));   %opts.lambda_val
    param_vals.lambda_b        = 2*opts.lambda_b;   % 4*;
    param_vals.lambda_history  = 0.09*opts.lambda_history; % 2*;
    param_vals.lambda_historyb = opts.lambda_historyb; %
    
    % Perform Recovery
    [A{ll}, ~] = BPDN_DF_bilinear(Z{ll}, MEAS_FUN, F, DWTfunc, param_vals);
%     [A2{ll}, ~] = BPDN_DF_bilinear(Z{ll}, MEAS_FUN, F_1, DWTfunc, param_vals);
    [A2{ll}, ~] = BPDN_DF_largescale(Z{ll}, MEAS_FUN, f_dyn, DWTfunc, param_vals);
    err_d(ll, :) = sum((X_ex{ll} - D*reshape(A{ll}{1}, size(D,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);
    err_nd(ll, :) = sum((X_ex{ll} - D*reshape(A2{ll}, size(D,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);
    fprintf('Sweep %d, Mean BPDF: %f, Mean D-BPDF: %f, Mean improvement: %f \n',ll, mean(mean(err_nd(1:ll, start_ss:end))), ...
        mean(mean(err_d(1:ll, start_ss:end))), mean(mean((err_nd(1:ll, start_ss:end) - err_d(1:ll, start_ss:end))./err_nd(1:ll, start_ss:end))))
% end

%%
e1 = reshape(err_d(1:ll, :), [], 1);
e2 = reshape(err_nd(1:ll, :), [], 1);

y_lims = [min([min(e1), min(e2)]), max([max(e1), max(e2)])];

nbins = 20;
bin_edges = logspace(log10(y_lims(1)), log10(y_lims(2)), nbins+1);

if floor(ll/1)*1 == ll
% figure(1)
% subplot(1,2,1), hist(log(reshape(err_d(1:ll, :), [], 1)), 50)
% title('Errors via Learned Dynamics')
% % set(gca, 'Xscale', 'log')
% subplot(1,2,2), hist(log(reshape(err_nd(1:ll, :), [], 1)), 50)
% title('Errors via BPDN-DF')
% % set(gca, 'Xscale', 'log')
% 
% figure(2)
% % subplot(2,1,1), plot([median(err_d, 1); mean(err_d, 1); median(err_nd, 1); mean(err_nd, 1)].')
% subplot(2,1,1), plot([mean(err_d(1:ll, :), 1); mean(err_nd(1:ll, :), 1)].')
% ylabel('Mean rMSE')
% xlabel('Time Index')
% legend('Learned Dynamics', 'BPDN-DF') 
% subplot(2,1,2), plot(mean((err_d(1:ll, :) - err_nd(1:ll, :))./err_nd(1:ll, :)))
% ylabel('rMSE % Improvement')
% xlabel('Time Index')

figure(3)
subplot(1,2,1), histogram(reshape(err_d(1:ll, :), [], 1), bin_edges,  'FaceAlpha', 0.3)
% hold on
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250])
title('Errors via Learned Dynamics', 'FontSize', 20)
subplot(1,2,2), histogram(reshape(err_nd(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250])
xlabel('rMSE', 'FontSize', 18)
% legend('Learned Dynamics', 'BPDN-DF')
% hold off
drawnow
else
end
end
% plot([sum((X_ex{x_sel} - D*reshape(A{1}, size(D,2), sig_opts.T_s)).^2,1)./sum(X_ex{x_sel}.^2);
%     sum((X_ex{x_sel} - D*reshape(A2{1}, size(D,2), sig_opts.T_s)).^2,1)./sum(X_ex{x_sel}.^2)]')
%%


e1 = reshape(err_d(1:ll, :), [], 1);
e2 = reshape(err_nd(1:ll, :), [], 1);

e1l = log(reshape(err_d(1:ll, :), [], 1));
e2l = log(reshape(err_nd(1:ll, :), [], 1));

y_lims = [min([min(e1), min(e2)]), max([max(e1), max(e2)])];

nbins = 20;
bin_edges = logspace(log10(y_lims(1)), log10(y_lims(2)), nbins+1);


figure(1)
subplot(1,2,1), histogram((reshape(err_d(1:ll, :), [], 1)), bin_edges)
title('Errors via Learned Dynamics')
set(gca, 'Xscale', 'log')
subplot(1,2,2), histogram((reshape(err_nd(1:ll, :), [], 1)), bin_edges)
title('Errors via BPDN-DF')
set(gca, 'Xscale', 'log')

figure(2)
% subplot(2,1,1), plot([median(err_d, 1); mean(err_d, 1); median(err_nd, 1); mean(err_nd, 1)].')
subplot(2,1,1), plot([mean(err_d(1:ll, :), 1); mean(err_nd(1:ll, :), 1)].')
ylabel('Mean rMSE')
xlabel('Time Index')
legend('Learned Dynamics', 'BPDN-DF') 
subplot(2,1,2), plot(mean((err_d(1:ll, :) - err_nd(1:ll, :))./err_nd(1:ll, :)))
ylabel('rMSE % Improvement')
xlabel('Time Index')


%%
figure(3)
subplot(1,2,1), histogram(reshape(err_d(1:ll, :), [], 1), bin_edges,  'FaceAlpha', 0.3)
% hold on
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600])
title('Errors via Learned Dynamics', 'FontSize', 20)
subplot(1,2,2), histogram(reshape(err_nd(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600])
xlabel('rMSE', 'FontSize', 18)
% legend('Learned Dynamics', 'BPDN-DF')
% hold off

%%

figure(3)
subplot(3,4,[1,10]), histogram(reshape(err_d(1:ll, :), [], 1), bin_edges,  'FaceAlpha', 0.3)
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600], 'Xtick', [0.01, 0.1, 1])
title('Learned Dynamics Errors', 'FontSize', 20)
line([mean(err_d(:)), mean(err_d(:))], [0,600], 'LineWidth', 3, 'LineStyle', '-', 'Color', [0,0.5,0])
line([median(err_d(:)), median(err_d(:))], [0,600], 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.5,0,0.5])
subplot(3,4,[3,12]), histogram(reshape(err_nd(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('BPDN-DF Errors', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600], 'Xtick', [0.01, 0.1, 1])
xlabel('rMSE', 'FontSize', 18)
line([mean(err_nd(:)), mean(err_nd(:))], [0,600], 'LineWidth', 3, 'LineStyle', '-', 'Color', [0,0.5,0])
line([median(err_nd(:)), median(err_nd(:))], [0,600], 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.5,0,0.5])
% subplot(3,4,[9,10,11,12]), plot([mean(err_d(1:ll, :), 1); mean(err_nd(1:ll, :), 1)].', 'LineWidth', 3)
% ylabel('Mean rMSE', 'FontSize', 18)
% xlabel('Time Index', 'FontSize', 18)
% legend('Learned Dynamics', 'BPDN-DF') 
% set(gca, 'FontSize', 16, 'Xlim', [1,20], 'Ylim', [0.14,0.27])


%%

figure(3)
subplot(3,4,[1,10]), histogram(reshape(err_d(1:ll, :), [], 1), 50,  'FaceAlpha', 0.3)
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600], 'Xtick', [0.01, 0.1, 1])
title('Learned Dynamics Errors', 'FontSize', 20)
line([mean(err_d(:)), mean(err_d(:))], [0,600], 'LineWidth', 3, 'LineStyle', '-', 'Color', [0,0.5,0])
line([median(err_d(:)), median(err_d(:))], [0,600], 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.5,0,0.5])
subplot(3,4,[3,12]), histogram(reshape(err_nd(1:ll, :), [], 1), 50, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('BPDN-DF Errors', 'FontSize', 20)
set(gca, 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600], 'Xtick', [0.01, 0.1, 1])
xlabel('rMSE', 'FontSize', 18)
line([mean(err_nd(:)), mean(err_nd(:))], [0,600], 'LineWidth', 3, 'LineStyle', '-', 'Color', [0,0.5,0])
line([median(err_nd(:)), median(err_nd(:))], [0,600], 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.5,0,0.5])