%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data
% use 10 latent states, 10 "recorded channels", 3000 timepoints
simulatedDmatrix = eye(10,10);
[latentStatesX, Fgt] = generateMultiSubNetwork(3000,5,5);
dFF = (simulatedDmatrix * latentStatesX).' + 0.0001*randn(size(latentStatesX,2), size(simulatedDmatrix,1));


%%

%dFF(dFF<0)=0; % remove negative values
dFF = dFF./max(abs(dFF),[],'all'); %set to abs
%%
figure(), cla;                                % ASC added 7/25
subplot(2,2,[1,2]), imagesc(dFF.')        % ASC added 7/25
subplot(2,2,3), imagesc(Fgt{1})           % ASC added 7/25
subplot(2,2,4), imagesc(Fgt{2})           % ASC added 7/25
%%

% Set parameters
inf_opts.nF              = 4; % 4
inf_opts.N               = 10; %scaled down from 100 (zebrafish, whole brain)
inf_opts.M               = size(dFF,2) ;
disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = 0.0001;   % ASC added 7/25 % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = 0.5;     % ASC added 7/25
inf_opts.lambda_b        = 0.025 ;       % ASC added 7/25 %0.025
inf_opts.lambda_historyb = 0.7;       % ASC added 7/25 %0.7
inf_opts.tol             = 1e-3; % 1e-3
inf_opts.max_iter2       = 20; %500 
inf_opts.max_iters       = 10;
inf_opts.special         = '';
inf_opts.F_update        = true;
inf_opts.D_update        = false;
inf_opts.N_ex            = 100; % 40
inf_opts.T_s             = 100; % 30
% inf_opts.special         = 'noobs';

inf_opts.step_d          = 20;  % ASC added 7/25
inf_opts.step_f          = 30;  % ASC added 7/25 % 30
inf_opts.plot_option     = 20;  % ASC added 7/25
inf_opts.lambda_f        = 0.05;

%%
% diary MSNDiary_3000i_10i2_Finitsvddecay_Kp5_Lval_p1_Lb_p1_tol_1em3_D10

[Phi, F] = bpdndf_dynamics_learning(dFF.', [], eye(10), inf_opts);              % Run the dictionary learning algorithm

% diary off

sizedFF = size(dFF.');
dFFcell = mat2cell(dFF.', sizedFF(1), sizedFF(2));
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients

subplot(2,1,1), plot(A_cell{1}.')   % ASC added 7/25
subplot(2,1,2), plot(B_cell{1}.')   % ASC added 7/25

% save('saveMSNWorkspace_3000i_10i2_Finitsvddecay_Kp5_Lval_p1_Lb_p1_tol_1em3_D10.mat'); 
%%

nTimePoints = size(latentStatesX,2);
figure();
subplot(3,1,1);
xNormLatent = zeros(nTimePoints,1);
for i = 1:nTimePoints
    xNormLatent(i) = norm(latentStatesX(:,i));
end
plot(xNormLatent);
title('Norm(Latent X)');
% colorbar;

subplot(3,1,2);
dFFNorm = zeros(nTimePoints,1);
dFF_t = dFF.';
for i = 1:nTimePoints
    dFFNorm(i) = norm(dFF_t(:,i));
end
plot(dFFNorm);
title('Norm(dFF)');
% colorbar;

subplot(3,1,3);
xNorm = zeros(nTimePoints,1);
for i = 1:nTimePoints
    xNorm(i) = norm(A_cell{1}(:,i));
end
plot(xNorm);
title('Norm(X)');
% colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
