%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data
% use 'sizeD' latent states, 'sizeD' "recorded channels", 'sizeT' timepoints

sizeD            = 8;
sizeT            = 3000; % T/2
simulatedDmatrix = eye(sizeD,sizeD);
contOpt          = false;
nF               = [3,3];

[latentStatesX, Fgt, groundTruthStates] = generateContMultiSubNetwork(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true);
dFF = statesToObservations(simulatedDmatrix, latentStatesX, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make ground truth

Fgt2 = cell(numel(Fgt{1})+numel(Fgt{2})+1,1);
for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(sizeD,sizeD);              end
for ll = 1:numel(Fgt{1}); Fgt2{ll}(1:sizeD/2,1:sizeD/2) = Fgt{1}{ll}; end
for ll = 1:numel(Fgt{2})
    Fgt2{ll+numel(Fgt{1})}(sizeD/2+1:sizeD,sizeD/2+1:sizeD) = Fgt{2}{ll}; 
end
Fgt2{end}                        = 0.001*randn(sizeD,sizeD);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Set parameters
inf_opts.nF              = 15; % # dynamics operators
inf_opts.N               = size(latentStatesX{1},1); % # latent states
inf_opts.M               = size(dFF{1},1) ; % original dimension (# channels)
% disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = 0.0001;         % ASC added 7/25 % VARY - tradeoff: approx. SNR - see ratios of traces % 0.0001 (don't shrink this further - gets to be essentially 0 for solver_L1RLS)
inf_opts.lambda_history  = 0.0001;         % ASC added 7/25
inf_opts.lambda_b        = 0.5;           % ASC added 7/25 %0.025
inf_opts.lambda_historyb = 0.45;           % ASC added 7/25 %0.7
inf_opts.tol             = 1e-8;           % 1e-3
inf_opts.max_iter2       = 20;             % 500 %20
inf_opts.max_iters       = 1000;
inf_opts.F_update        = true; % default = true;
inf_opts.D_update        = false; % NoObs case - default = true;
inf_opts.N_ex            = 50;   
inf_opts.T_s             = 200; 
inf_opts.step_d          = 1;  % ASC added 7/25
inf_opts.step_f          = 10;  % ASC added 7/25 % 30
inf_opts.step_decay      = 0.998;  
inf_opts.plot_option     = 10;  % ASC added 7/25
inf_opts.lambda_f_decay  = 0.996;
inf_opts.lambda_f        = 0.20;
inf_opts.solver_type     = ''; % fista
inf_opts.special         = 'noobs';
inf_opts.deltaDynamics   = false; %default: false %use x_t-x_{t-1}

% inf_opts.solver_type     = 'tfocs'; %default: ''; (fista)
% inf_opts.CVX_Precision   = 'default';   %default = 'default'
% inf_opts.special         = ''; % regular bilinear inference
% inf_opts.debias          = true; % default = true; not used for NoObs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Phi, F] = bpdndf_dynamics_learning(dFF, [], simulatedDmatrix, inf_opts);              % Run the dictionary learning algorithm
% Phi = simulatedDmatrix; F = Fgt2;
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients

%

plotMultiNetworkOutput(5, dFF, A_cell, B_cell, groundTruthStates, F)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%