%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RWL1-DF Testing script
%
% This test script loads and subsamples a number of different signal 
% types, including video, MRI and communications channel data. The 
% resulting samples are used to recover the sequence using basis pursuit 
% de-noising (BPDN), re-weighted \ell_1 (RWL1), BPDN dynamic filtering 
% (BPDN-DF) and RWL1 dynamic filtering (RWL1-DF).
% 
% 
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated August 20, 2012. 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Parameters

% Running Options 
reset_samp = 1;
same_samp = 0;
calc_l1 = 1;
calc_rwl1 = 1;
calc_dl1 = 1;
calc_drwl1 = 1;
calc_dcsamp = 0;
calc_modcs = 1;
calc_rwl2 = 0;
calc_wl1p = 1;
calc_drwl1b = 1;
normalize_opt = 0;

num_algs = 8;

sim_type = 'mri';

sim_param_setup

% General Parameters
TOL = 1e-3;          % TFOCS tolerance parameter

% Initialize result arrays
PSNR_ALL = zeros(load_params.T_s, num_algs, num_trials);
rMSE_ALL = zeros(load_params.T_s, num_algs, num_trials);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data

for trial_index = 1:num_trials

    if reset_samp == 1
        [y,MEAS_FUN,x0] = Load_Data(sim_type,filename, load_params);
    end    
	    
    %% BPDN Reconstruction
    if calc_l1 == 1
        lambda_val = 0.0125;
        [coef_cs, recon_cs, rMSE_cs, PSNR_cs] = ...
               BPDN_largescale(y, MEAS_FUN, lambda_val, TOL, DWTfunc, x0);
    end

    %% RWL1 Reconstruction
    if calc_rwl1 == 1
        % 0.3400    0.1200    0.0016
        lambda_val = 0.0016; % 0.001;
        rwl1_reg = 0.12;     % 0.1;
        rwl1_mult = 0.34;    % 0.05;
        [coef_rwcs, recon_rwcs, rMSE_rwcs, PSNR_rwcs, conv_vals_store2] = ...
              RWL1_largescale(y, MEAS_FUN, [lambda_val, rwl1_reg, rwl1_mult], ...
              TOL, DWTfunc, x0);
    end

    %% BPDN-DF Reconstruction
    if calc_dl1 == 1
        param_vals.lambda_val = 0.017;   % 0.01;
        param_vals.lambda_history = 0.01; % 0.3;
        param_vals.tol = TOL;
        [coef_dcs, recon_dcs, rMSE_dcs, PSNR_dcs] = ...
            BPDN_DF_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);
    end
    
    %% RWL1-DF Reconstruction
    if calc_drwl1 == 1
        %  0.1200    0.2300    0.0129 
        param_vals.lambda_val = 0.0129;  % 0.001;
        param_vals.rwl1_reg = 0.23;      % 0.2;
        param_vals.rwl1_mult = 0.12;     % 0.4;
        param_vals.beta = 1;
        param_vals.tol = TOL;
        [coef_drw, recon_drw, rMSE_drw, PSNR_drw, conv_vals_store] = ...
            RWL1_DF_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);
    end
    
    %% DCS-AMP Reconstruction
    if calc_dcsamp == 1
%         DCS_params.lambda_0 = 0.1;   % 0.17
%         DCS_params.pz1      = 0.04;   % 0.03
%         DCS_params.p1z      = 0.1;    % 0.1
%         DCS_params.eta_0    = 0.01;      % 0
%         DCS_params.kappa_0  = 1;      % 1
%         DCS_params.alpha    = 0.001;  % 0.001    
%         DCS_params.rho      = 0.2;    % 0.1, 1e-5;      
%         DCS_params.sig2e    = 0.002;  % 0.002  
%         DCS_params.eps      = 1e-5;   % 1e-5
%         
        DCS_params.lambda_0 = 0.10;     % 0.10
        DCS_params.pz1      = 0.04;     % 0.03
        DCS_params.p1z      = DCS_params.pz1 * DCS_params.lambda_0 / (1 - DCS_params.pz1);   % Steady-state sparsity of lambda_0
        DCS_params.eta_0    = 0;        % 0
        DCS_params.kappa_0  = 1;        % 1
        DCS_params.alpha    = 0.10;     % 0.001    
        DCS_params.rho      = DCS_params.kappa_0 / DCS_params.alpha^2;	% Steady-state variance of kappa_0    
        DCS_params.sig2e    = 0.001;    % 0.002  
        DCS_params.eps      = 1e-5;     % 1e-5

        DCS_options.eq_iter = 25;       % Inner AMP iterations per turbo iteration
        DCS_options.verbose = true;     % Print run-time information
	    
        % causal without normalization
        DCS_options.smooth_iter = -1;
        DCS_options.update = 0;
%         [vid_coef_amp, vid_recon_amp, vid_rMSE_amp, vid_PSNR_amp] = ...
% 	    DCSAMP_largescale(y, MEAS_FUN, DWTfunc, DCS_params, [], DCS_options, x0);
        % causal with normalization
        [coef_amp2, recon_amp2, rMSE_amp2, PSNR_amp2] = ...
	    DCSAMP_largescale(y, MEAS_FUN, DWTfunc, DCS_params, NORMS_VEC, DCS_options, x0);
    end

    %% modCS Reconstruction
    if calc_modcs == 1
        param_vals.lambda_val = 0.3816;
        param_vals.EPS = 0.0013*sum(y(:, :, 1).^2);
        param_vals.tol = TOL;
        [coef_modcs, recon_modcs, rMSE_modcs, PSNR_modcs] = ...
            modCS_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);
    end

    %% RWL2 Reconstruction
    if calc_rwl2 == 1
        param_vals.lambda_val = 0.001;
        param_vals.rwl1_reg = 0.01;
        %param_vals.rwl1_mult = 0.4; %0.4
        param_vals.beta = 10; %0.01;
        param_vals.tol = TOL;
        [coef_rwl2b, recon_rwl2b, rMSE_rwl2b, PSNR_rwl2b] = ...
            RWL2mod_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);
    end
    
    %% RWL1-DF Multiplication Reconstruction
    if calc_drwl1b == 1
        param_vals.lambda_val = 0.001;
        param_vals.rwl1_reg = 0.24;
        param_vals.rwl1_mult = 1;
        param_vals.beta = 0.4;
        param_vals.tol = TOL;
        [coef_mdrw, recon_mdrw, rMSE_mdrw, PSNR_mdrw] = ...
            RWL1_DF2_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);
    end
 
    if calc_wl1p == 1
        param_vals.lambda_val = 0.01;
        param_vals.rwl1_reg = 1;
        param_vals.rwl1_mult = 1;
        param_vals.beta = 0.55;
        param_vals.tol = TOL;
        [coef_wl1p, recon_wl1p, rMSE_wl1p, PSNR_wl1p] = ...
            WL1_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);
    end

   
    %% Aggregate the results:
    PSNR_ALL(:, :, trial_index) = [PSNR_cs(:), PSNR_rwcs(:), PSNR_dcs(:), PSNR_drw(:), PSNR_amp2(:), PSNR_mdrw(:), PSNR_modcs(:), PSNR_wl1p(:)];
    rMSE_ALL(:, :, trial_index) = [rMSE_cs(:), rMSE_rwcs(:), rMSE_dcs(:), rMSE_drw(:), rMSE_amp2(:), rMSE_mdrw(:), rMSE_modcs(:), rMSE_wl1p(:)];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting

% save DTDWT_wCom.mat PSNR_ALL rMSE_ALL M N

alg_names = {'BPDN', 'RWL1', 'BPDN-DF', 'RWL1-DF'};
bin_cents = linspace(min(reshape(rMSE_ALL, [], 1)), max(reshape(rMSE_ALL, [], 1)), 30);

% Time plot
line_style_list = {':b', '-g', '-.r', '--c', ':k', '--m', '-y', '--m'};
line_style_list2 = {'.b', '.g', '.r', '.c', '.k', '^m', 'oy', '^m'};
color_style_list = {'b', 'g', 'r', 'c', 'k', 'm', 'y', 'm'};

figure, hold on;
for kk = 1:8
    plot(mean(rMSE_ALL(:, kk), 3), line_style_list{kk}, 'LineWidth', 3)
%     plot(rMSE_ALL(:, kk, 2), line_style_list{kk}, 'LineWidth', 3)
end
box on
legend('Independent BPDN', 'Independent rw-BPDN', 'Dynamic BPDN', 'Dynamic rw-BPDN', 'DCS-AMP', 'RWL1-SS', 'modCS', 'WL1P')
set(gca, 'FontSize', 18, 'Xlim', [1,size(rMSE_ALL, 1)], 'Ylim', [0, 0.1])
xlabel('Frame Number', 'FontSize', 22)
ylabel('rMSE', 'FontSize', 22)

%% Histogram plot
FULL_MEANS = mean(mean(rMSE_ALL, 3), 1);
FULL_MEDIANS = median(mean(rMSE_ALL, 3), 1);
FIG_HEIGHT = 150;
Y_LIMS_PER1 = [35, 55, 35,60];
alg_names = {'BPDN', 'RWL1', 'BPDN-DF', 'RWL1-DF'};
bin_cents = linspace(min(min(rMSE_ALL)), max(max(rMSE_ALL)), 30);
figure;
for kk = 1:4
    subplot(2, 2, kk), hold on
    hist(mean(rMSE_ALL(:, kk, :), 3), bin_cents)
    line([FULL_MEANS(kk), FULL_MEANS(kk)], [0, FIG_HEIGHT], 'LineStyle', ...
        '--', 'LineWidth', 3, 'Color', [0,0,0])
    box on
    title(alg_names{kk}, 'FontSize', 24)
    set(gca, 'FontSize', 20, 'YLim', [0, Y_LIMS_PER1(kk)], ...
        'XLim', [0, 0.065], 'XTick', [0.02,0.04,0.06])
    xlabel('rMSE', 'FontSize', 24)
    axis square
    plotarrow('Start', [FULL_MEDIANS(kk), -15*Y_LIMS_PER1(kk)/150], 'Stop', ...
        [FULL_MEDIANS(kk), 0], 'Length', 4, 'BaseAngle', 85, ...
        'TipAngle', 35, 'Width', 3)
    hold off
end

% %% DCS AMP test plots
% 
% DCS_rMSE_ALL = [vid_rMSE_amp(:), vid_rMSE_amp2(:), vid_rMSE_amp3(:), vid_rMSE_amp4(:), vid_rMSE_amp5(:)];
% 
% alg_names = {'Norm', 'No norm', 'Smooth update w/ norm', 'Smooth update w/o norm', 'Smooth no update w/ norm'};
% bin_cents = linspace(min(min(rMSE_ALL)), max(max(rMSE_ALL)), 30);
% 
% % Time plot
% line_style_list = {':b', '-g', '-.r', '--c', ':k'};
% line_style_list2 = {'.b', '.g', '.r', '.c', '.k'};
% color_style_list = {'b', 'g', 'r', 'c', 'k'};
% 
% figure, hold on;
% for kk = 1:5
%     plot(mean(DCS_rMSE_ALL(:, kk), 3), line_style_list{kk}, 'LineWidth', 3)
% end
% box on
% legend('No Norm', 'Norm', 'Smooth update w/ norm', 'Smooth update w/o norm', 'Smooth no update w/ norm')
% set(gca, 'FontSize', 18, 'Xlim', [1,size(rMSE_ALL, 1)], 'Ylim', [0.012, 0.095])
% xlabel('Frame Number', 'FontSize', 22)
% ylabel('Mean rMSE', 'FontSize', 22)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
