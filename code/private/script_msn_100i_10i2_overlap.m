1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% [latentStatesX, Fgt] = generateMultiSubNetwork_Overlap(3000,5,5);
[latentStatesX] = generateMultiSubNetwork_Overlap(3000,5,5);
dFF = (simulatedDmatrix * latentStatesX).' + 0.0001*randn(size(latentStatesX,2), size(simulatedDmatrix,1));


%%

%dFF(dFF<0)=0; % remove negative values
dFF = dFF./max(abs(dFF),[],'all'); %set to abs

figure(1);                                % ASC added 7/25
subplot(2,2,[1,2]), imagesc(dFF.')        % ASC added 7/25
% subplot(2,2,3), imagesc(Fgt{1})           % ASC added 7/25
% subplot(2,2,4), imagesc(Fgt{2})           % ASC added 7/25
%%

% Set parameters
inf_opts.nF              = 4;
inf_opts.N               = 10; %scaled down from 100 (zebrafish, whole brain)
inf_opts.M               = size(dFF,2) ;
disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = 0.0001;   % ASC added 7/25 % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = 0.5;     % ASC added 7/25
inf_opts.lambda_b        = 0.025;       % ASC added 7/25
inf_opts.lambda_historyb = 0.7;       % ASC added 7/25
inf_opts.tol             = 1e-3; % 1e-3
inf_opts.max_iter2       = 20; %500 
inf_opts.max_iters       = 20;
inf_opts.special         = '';
inf_opts.F_update        = true;
inf_opts.D_update        = false;
inf_opts.N_ex            = 40;
inf_opts.T_s             = 30;
% inf_opts.special         = 'noobs';

inf_opts.step_d          = 20;  % ASC added 7/25
inf_opts.step_f          = 30;  % ASC added 7/25
inf_opts.plot_option     = 20;  % ASC added 7/25
inf_opts.lambda_f        = 0.05;

%%
diary MSNDiary_100i_10i2_overlap
[Phi, F] = bpdndf_dynamics_learning(dFF.', [], eye(10), inf_opts);              % Run the dictionary learning algorithm

diary off

sizedFF = size(dFF.');
dFFcell = mat2cell(dFF.', sizedFF(1), sizedFF(2));
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients

subplot(2,1,1), plot(A_cell{1}.')   % ASC added 7/25
subplot(2,1,2), plot(B_cell{1}.')   % ASC added 7/25

save('saveMSNWorkspace_100i_10i2_overlap.mat');               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
