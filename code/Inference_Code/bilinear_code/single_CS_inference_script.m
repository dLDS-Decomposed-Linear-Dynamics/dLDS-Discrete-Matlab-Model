%% Mini Inference script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data 
clear

load /home/adam/GITrepos/dynamics_learning/results/20180521_BBCsysid_all.mat
load /home/adam/GITrepos/dynamics_learning/results/OLD/20150216_BBCvideo_12x12_4x_nF20.mat 
% load /home/adam/GITrepos/dynamics_learning/results/OLD/20150218_BBCvideo_12x12_4x_nF20.mat
F_20 = F;
D_20 = D;
D_nd = D_bbcsid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample data

x_sel  = 20;
n_meas = ceil(0.2*size(D,1));
n_var = 0.001; % 0.0005;

X_ex = cell(1, x_sel);                                                     % Initialize the cell array of videos to test
for ll = 1:x_sel                                                           % Get a random video
    while sum(sum(abs(X_ex{ll}))==0)>0
        X_ex{ll} = rand_bbc_video([sqrt(sig_opts.M), sqrt(sig_opts.M),...
                                           sig_opts.T_s], sig_opts.dsamp); % Sample a few video sequences
%   X_ex{kk} = rand_bbc_video([sqrt(inf_opts.M), ... 
%                          sqrt(inf_opts.M), inf_opts.T_s], inf_opts.dsamp); % Sample some videos from the BBC dataset
    end
end

G = cell(1,numel(X_ex));                                                   % Initialize cell array of sensing matrices
Z = cell(1,numel(X_ex));                                                   % Initialize cell array of measurements (X_ex passed through the sensing matrix)

for ll = 1:numel(X_ex)
    for kk = 1:sig_opts.T_s
        Gtmp = randn(n_meas, size(D, 1));                                  % Make a random sensing matrix
        Gtmp = Gtmp*diag(1./sqrt(sum(Gtmp.^2,1)));                         % Normalize sensing matrix appropriately
        G{ll}(:,:,kk) = Gtmp;                                              % Store the current sensing matrix
        Z{ll}(:,:,kk) = G{ll}(:,:,kk)*X_ex{ll}(:,kk)+n_var*randn(n_meas,1);% Create noisy samples
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

opts                         = inf_opts;
param_vals_nd.lambda_val     = 0.05*opts.lambda_val;      % 8*;
param_vals_nd.lambda_history = 0.5*opts.lambda_history;       % 2*;

param_vals_si.lambda_val      = 0.01*opts.lambda_val;      % 8*;
param_vals_si.lambda_history  = 2*opts.lambda_history;       % 2*;

param_vals_1.lambda_val      = 0.1*opts.lambda_val;      % 8*;
param_vals_1.lambda_b        = opts.lambda_b;             % 4*
param_vals_1.lambda_history  = opts.lambda_history;       % 2*;

param_vals_6.lambda_val      = 0.01*opts.lambda_val;      % 8*;
param_vals_6.lambda_b        = opts.lambda_b;             % 4*
param_vals_6.lambda_history  = opts.lambda_history;       % 2*;

param_vals_12.lambda_val      = 0.01*opts.lambda_val;      % 8*;
param_vals_12.lambda_b        = opts.lambda_b;             % 4*
param_vals_12.lambda_history  = opts.lambda_history;       % 2*;

param_vals_20.lambda_val      = 0.01*opts.lambda_val;      % 8*;
param_vals_20.lambda_b        = 4*opts.lambda_b;             % 4*
param_vals_20.lambda_history  = 2*opts.lambda_history;       % 2*;

if isfield(opts,'tol')
    param_vals.tol = opts.tol;
else
    param_vals.tol = 1e-3;
end

MEAS_FUN_nd    = cell(1,sig_opts.T_s);                                     % Initialize the arrays of measurement functions for all methods
MEAS_FUN_si    = cell(1,sig_opts.T_s);                                     %   |
MEAS_FUN_1     = cell(1,sig_opts.T_s);                                     %   |
MEAS_FUN_6     = cell(1,sig_opts.T_s);                                     %   |
MEAS_FUN_12    = cell(1,sig_opts.T_s);                                     %   |
MEAS_FUN_20    = cell(1,sig_opts.T_s);                                     % -----
DWTfunc.apply  = @(q) q;                                                   % For now set the wavelet matrix to identity
DWTfunc.invert = @(q) q;                                                   % For now set the wavelet matrix to identity



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize error and coefficient arrays

A_nd = cell(1, numel(X_ex));
A_si = cell(1, numel(X_ex));
A_1  = cell(1, numel(X_ex));
A_6  = cell(1, numel(X_ex));
A_12 = cell(1, numel(X_ex));
A_20 = cell(1, numel(X_ex));

err_nd = zeros(numel(X_ex), sig_opts.T_s);
err_si = zeros(numel(X_ex), sig_opts.T_s);
err_1  = zeros(numel(X_ex), sig_opts.T_s);
err_6  = zeros(numel(X_ex), sig_opts.T_s);
err_12 = zeros(numel(X_ex), sig_opts.T_s);
err_20 = zeros(numel(X_ex), sig_opts.T_s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform inference

start_ss = 1;

for ll = 1:numel(X_ex)
    
    for kk = 1:sig_opts.T_s
        meas_func.Phi = @(z) G{ll}(:,:,kk)*D_nd*z;
        meas_func.Phit = @(z) ((G{ll}(:,:,kk)*D_nd)')*z;
        MEAS_FUN_nd{kk} = meas_func;

        meas_func.Phi = @(z) G{ll}(:,:,kk)*D_bbcsid*z;
        meas_func.Phit = @(z) ((G{ll}(:,:,kk)*D_bbcsid)')*z;
        MEAS_FUN_si{kk} = meas_func;
        
        meas_func.Phi = @(z) G{ll}(:,:,kk)*D_1*z;
        meas_func.Phit = @(z) ((G{ll}(:,:,kk)*D_1)')*z;
        MEAS_FUN_1{kk} = meas_func;

        meas_func.Phi = @(z) G{ll}(:,:,kk)*D_6*z;
        meas_func.Phit = @(z) ((G{ll}(:,:,kk)*D_6)')*z;
        MEAS_FUN_6{kk} = meas_func;

        meas_func.Phi = @(z) G{ll}(:,:,kk)*D_12*z;
        meas_func.Phit = @(z) ((G{ll}(:,:,kk)*D_12)')*z;
        MEAS_FUN_12{kk} = meas_func;

        meas_func.Phi = @(z) G{ll}(:,:,kk)*D_20*z;
        meas_func.Phit = @(z) ((G{ll}(:,:,kk)*D_20)')*z;
        MEAS_FUN_20{kk} = meas_func;
    end
    
%     % Modify Parameters?
%     param_vals.lambda_val      = 0.01*opts.lambda_val;  % max(max((G{ll}(:,:,1)')*Z{ll}(:,:,1).^2));   %opts.lambda_val
%     param_vals.lambda_b        = opts.lambda_b;         % 4*;
%     param_vals.lambda_history  = opts.lambda_history;   % 2*;
%     param_vals.lambda_historyb = opts.lambda_historyb;  %
    
    % Perform Recovery
    [A_nd{ll}, ~] = BPDN_DF_largescale(Z{ll}, MEAS_FUN_nd, @(q) q, DWTfunc, param_vals_nd);
    [A_si{ll}, ~] = BPDN_DF_largescale(Z{ll}, MEAS_FUN_si, @(q) F_bbcsid{1}*q, DWTfunc, param_vals_si);
    [A_1{ll}, ~]  = BPDN_DF_bilinear(Z{ll}, MEAS_FUN_1, F_1, DWTfunc, param_vals_1);
    [A_6{ll}, ~]  = BPDN_DF_bilinear(Z{ll}, MEAS_FUN_6, F_6, DWTfunc, param_vals_6);
    [A_12{ll}, ~] = BPDN_DF_bilinear(Z{ll}, MEAS_FUN_12, F_12, DWTfunc, param_vals_12);
    [A_20{ll}, ~] = BPDN_DF_bilinear(Z{ll}, MEAS_FUN_20, F_20, DWTfunc, param_vals_20);
    
    
    err_nd(ll, :) = sum((X_ex{ll} - D_nd*reshape(A_nd{ll}, size(D_nd,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);
    err_si(ll, :) = sum((X_ex{ll} - D_bbcsid*reshape(A_si{ll}, size(D_bbcsid,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);
    err_1(ll, :)  = sum((X_ex{ll} - D_1*reshape(A_1{ll}{1}, size(D_1,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);
    err_6(ll, :)  = sum((X_ex{ll} - D_6*reshape(A_6{ll}{1}, size(D_6,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);
    err_12(ll, :) = sum((X_ex{ll} - D_12*reshape(A_12{ll}{1}, size(D_12,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);
    err_20(ll, :) = sum((X_ex{ll} - D_20*reshape(A_20{ll}{1}, size(D_20,2), sig_opts.T_s)).^2,1)./sum(X_ex{ll}.^2);

    e_nd = reshape(err_nd(1:ll, :), [], 1);
    e_si = reshape(err_si(1:ll, :), [], 1);
    e_1  = reshape(err_1(1:ll, :), [], 1);
    e_6  = reshape(err_6(1:ll, :), [], 1);
    e_12 = reshape(err_12(1:ll, :), [], 1);
    e_20 = reshape(err_20(1:ll, :), [], 1);
    
    fprintf('Sweep %d: \n', ll)
    fprintf('\t\t BPDF: Mean = %f, median = %f \n', mean(e_nd), median(e_nd))
    fprintf('\t\t sysID: Mean = %f, median = %f, pct improvement = %f \n', mean(e_si), median(e_si), mean((e_nd - e_si)./e_nd))
    fprintf('\t\t L-BPDF(1): Mean = %f, median = %f, pct improvement = %f \n', mean(e_1), median(e_1), mean((e_nd - e_1)./e_nd))
    fprintf('\t\t L-BPDF(6): Mean = %f, median = %f, pct improvement = %f \n', mean(e_6), median(e_6), mean((e_nd - e_6)./e_nd))
    fprintf('\t\t L-BPDF(12): Mean = %f, median = %f, pct improvement = %f \n', mean(e_12), median(e_12), mean((e_nd - e_12)./e_nd))
    fprintf('\t\t L-BPDF(25): Mean = %f, median = %f, pct improvement = %f \n',mean(e_20), median(e_20), mean((e_nd - e_20)./e_nd))
% end

    %



    y_lims    = [min([min(e_nd), min(e_si), min(e_1), min(e_6), min(e_12), min(e_20)]), ...
                      max([max(e_nd), max(e_si), max(e_1), max(e_6), max(e_12), max(e_20)])];
    nbins     = 20;
    bin_edges = logspace(log10(y_lims(1)), log10(y_lims(2)), nbins+1);

    if mod(ll,1) == 0
        figure(1)
        subplot(2,3,1), histogram(reshape(err_nd(1:ll, :), [], 1), bin_edges,  'FaceAlpha', 0.3)
        xlabel('rMSE', 'FontSize', 18)
        set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
        title('Errors via BPDN-DF', 'FontSize', 20)

        subplot(2,3,2), histogram(reshape(err_si(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
        title('Errors via BPDN-DF', 'FontSize', 20)
        set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
        xlabel('rMSE', 'FontSize', 18)
        title('Errors via sysID Dynamics', 'FontSize', 20)

        subplot(2,3,3), histogram(reshape(err_1(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
        title('Errors via BPDN-DF', 'FontSize', 20)
        set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
        xlabel('rMSE', 'FontSize', 18)
        title('Errors via 1 Learned Dynamics', 'FontSize', 20)

        subplot(2,3,4), histogram(reshape(err_6(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
        title('Errors via BPDN-DF', 'FontSize', 20)
        set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
        xlabel('rMSE', 'FontSize', 18)
        title('Errors via 6 Learned Dynamics', 'FontSize', 20)

        subplot(2,3,5), histogram(reshape(err_12(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
        title('Errors via BPDN-DF', 'FontSize', 20)
        set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
        xlabel('rMSE', 'FontSize', 18)
        title('Errors via 12 Learned Dynamics', 'FontSize', 20)

        subplot(2,3,6), histogram(reshape(err_20(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
        title('Errors via BPDN-DF', 'FontSize', 20)
        set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
        xlabel('rMSE', 'FontSize', 18)
        title('Errors via 25 Learned Dynamics', 'FontSize', 20)
        set(gcf,'color',[1,1,1])
        drawnow
    else
    end
end


%%


e_nd = reshape(err_nd(1:ll, :), [], 1);
e_si = reshape(err_si(1:ll, :), [], 1);
e_1  = reshape(err_1(1:ll, :), [], 1);
e_6  = reshape(err_6(1:ll, :), [], 1);
e_12 = reshape(err_12(1:ll, :), [], 1);
e_20 = reshape(err_20(1:ll, :), [], 1);

el_nd = log(reshape(err_nd(1:ll, :), [], 1));
el_si = log(reshape(err_si(1:ll, :), [], 1));
el_1  = log(reshape(err_1(1:ll, :), [], 1));
el_6  = log(reshape(err_6(1:ll, :), [], 1));
el_12 = log(reshape(err_12(1:ll, :), [], 1));
el_25 = log(reshape(err_20(1:ll, :), [], 1));

y_lims = [min([min(e_nd), min(e_si), min(e_1), min(e_6), min(e_12), min(e_20)]), ...
                  max([max(e_nd), max(e_si), max(e_1), max(e_6), max(e_12), max(e_20)])];

nbins = 20;
bin_edges = logspace(log10(y_lims(1)), log10(y_lims(2)), nbins+1);


figure(1)
subplot(2,3,1), histogram(reshape(err_nd(1:ll, :), [], 1), bin_edges,  'FaceAlpha', 0.3)
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
title('Errors via BPDN-DF', 'FontSize', 20)

subplot(2,3,2), histogram(reshape(err_si(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
xlabel('rMSE', 'FontSize', 18)
title('Errors via sysID Dynamics', 'FontSize', 20)

subplot(2,3,3), histogram(reshape(err_1(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
xlabel('rMSE', 'FontSize', 18)
title('Errors via 1 Learned Dynamics', 'FontSize', 20)

subplot(2,3,4), histogram(reshape(err_6(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
xlabel('rMSE', 'FontSize', 18)
title('Errors via 6 Learned Dynamics', 'FontSize', 20)

subplot(2,3,5), histogram(reshape(err_12(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
xlabel('rMSE', 'FontSize', 18)
title('Errors via 12 Learned Dynamics', 'FontSize', 20)

subplot(2,3,6), histogram(reshape(err_20(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,250],'XTick',[0.01,0.1,1])
xlabel('rMSE', 'FontSize', 18)
title('Errors via 25 Learned Dynamics', 'FontSize', 20)
set(gcf,'color',[1,1,1])
drawnow

figure(2)
% subplot(2,1,1), plot([median(err_d, 1); mean(err_d, 1); median(err_nd, 1); mean(err_nd, 1)].')
subplot(2,1,1), plot([mean(err_nd(1:ll, :), 1); mean(err_si(1:ll, :), 1)].')
ylabel('Mean rMSE')
xlabel('Time Index')
legend('Learned Dynamics', 'BPDN-DF') 
subplot(2,1,2), plot(mean((err_d(1:ll, :) - err_nd(1:ll, :))./err_nd(1:ll, :)))
ylabel('rMSE % Improvement')
xlabel('Time Index')


%%
figure(3)
subplot(1,2,1), histogram(reshape(err_nd(1:ll, :), [], 1), bin_edges,  'FaceAlpha', 0.3)
% hold on
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600])
title('Errors via Learned Dynamics', 'FontSize', 20)
subplot(1,2,2), histogram(reshape(err_20(1:ll, :), [], 1), bin_edges, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('Errors via BPDN-DF', 'FontSize', 20)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600])
xlabel('rMSE', 'FontSize', 18)
% legend('Learned Dynamics', 'BPDN-DF')
% hold off

%%

figure(3)
subplot(3,4,[1,10]), histogram(reshape(err_20(1:ll, :), [], 1), bin_edges,  'FaceAlpha', 0.3)
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xscale', 'log', 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600], 'Xtick', [0.01, 0.1, 1])
title('Learned Dynamics Errors', 'FontSize', 20)
line([mean(err_20(:)), mean(err_20(:))], [0,600], 'LineWidth', 3, 'LineStyle', '-', 'Color', [0,0.5,0])
line([median(err_20(:)), median(err_20(:))], [0,600], 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.5,0,0.5])
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
subplot(3,4,[1,10]), histogram(reshape(err_20(1:ll, :), [], 1), 50,  'FaceAlpha', 0.3)
xlabel('rMSE', 'FontSize', 18)
set(gca, 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600], 'Xtick', [0.01, 0.1, 1])
title('Learned Dynamics Errors', 'FontSize', 20)
line([mean(err_20(:)), mean(err_20(:))], [0,600], 'LineWidth', 3, 'LineStyle', '-', 'Color', [0,0.5,0])
line([median(err_20(:)), median(err_20(:))], [0,600], 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.5,0,0.5])
subplot(3,4,[3,12]), histogram(reshape(err_nd(1:ll, :), [], 1), 50, 'FaceAlpha', 0.5, 'FaceColor', 'red')
title('BPDN-DF Errors', 'FontSize', 20)
set(gca, 'Xlim', [0.9*y_lims(1), 1.1*y_lims(2)], 'FontSize', 16, 'Ylim', [0,600], 'Xtick', [0.01, 0.1, 1])
xlabel('rMSE', 'FontSize', 18)
line([mean(err_nd(:)), mean(err_nd(:))], [0,600], 'LineWidth', 3, 'LineStyle', '-', 'Color', [0,0.5,0])
line([median(err_nd(:)), median(err_nd(:))], [0,600], 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.5,0,0.5])


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