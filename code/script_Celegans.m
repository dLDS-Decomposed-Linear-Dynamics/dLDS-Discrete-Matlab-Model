%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parpool(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data (as a cell array)

dat = load('./CElegansData/WT_Stim.mat');

dFF_S_AllWorms = {dat.WT_Stim.traces}.';
dFF_S_AllWorms = cellfun(@(x) x.',dFF_S_AllWorms,UniformOutput=false); % works in Matlab 2022a but not 2020a
bhv_S_AllWorms = {dat.WT_Stim.States};

%% create cell array with your samples 

%Examples

% % 1 worm, 1 trial
% whichStimWorm = 1;
% dFF = dFF_S_AllWorms(whichStimWorm);

% 1 worm, multiple durations of trials (artificially trimmed) (same channels, could be different trial durations)
whichStimWorm = 1;
dFF = {dFF_S_AllWorms{whichStimWorm}(:,1:2000);dFF_S_AllWorms{whichStimWorm}(:,1:1000);dFF_S_AllWorms{whichStimWorm}(:,1:500)}; 
bhv = {bhv_S_AllWorms{whichStimWorm}(:,1:2000);bhv_S_AllWorms{whichStimWorm}(:,1:1000);bhv_S_AllWorms{whichStimWorm}(:,1:500)}; 


%% Preprocess data
% Because of the expected distribution of calcium imaging data for 
% C. elegans, we replaced negative values with zeros. This affected up to
% 3 percent of the values for any particular worm.
% We then rescaled the data between 0 and 1 across the whole dataset, not
% neuron by neuron, because we wanted to preserve the relationships between
% neurons, with less active neurons being represented accordingly in the
% latent space.

for ii = 1:size(dFF,1)
    thisDFF = dFF{ii};
    thisDFF(thisDFF<0) = 0;
    thisDFF = thisDFF./max(thisDFF,[],'all');
    dFF{ii} = thisDFF;
end


%%

% Set parameters - must scale up to the data range
% In our case we
inf_opts.nF              = 10; % dynamics operators
inf_opts.N               = 10; % latent dimensions
% Closer to 0 means weaker regularization - that term becomes less
% important in the optimization.
inf_opts.lambda_val      = 0.1; % regularization: sparsity of states (x), i.e., ||x_t||
inf_opts.lambda_history  = 0.9; % regularization: smoothness of latent states (x) from time point to time point, i.e., ||x_t-sum^M(f_{m}c_{mt})x_{t-1}||
inf_opts.lambda_b        = 0.01; % regularization: sparsity of dynamics coefficients (c), i.e., ||c_t||
inf_opts.lambda_historyb = 0; % regularization: smoothness of dynamics coefficients (c), i.e., ||c-c_{t-1}||
inf_opts.tol             = 1e-3; 
inf_opts.max_iter2       = 50; 
inf_opts.max_iters       = 500; % parameter that controls number of EM iterations run
inf_opts.special         = ''; % default
inf_opts.F_update        = true; % true: update dynamics operators
inf_opts.D_update        = true; % if inf_opts.special is "noobs", set this to false (no observations, no observation matrix)
inf_opts.N_ex            = 40; % number of examples selected from the data
inf_opts.T_s             = 30; % number of time points in each example


inf_opts.M = size(dFF{1},1) ; % original dimension (# channels)


%%

[Phi, F] = bpdndf_dynamics_learning(dFF, [], [], inf_opts);              % Run the dictionary learning algorithm

[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients

%% Display variance explained overall for the reconstruction
stdzdOn = false; % true if you want to rescale the reconstructed data to the original range (not necessary)

for ii = 1:size(A_cell,1) %in the event of multiple trials
    x_values = A_cell{ii};
    dataReconstructed = (Phi * x_values).'; % from inferred latent states and learned dictionary    
    
    maxvalDFF = max(dFF{ii},[],'all');
    
    if stdzdOn 
        dataReconstructedRescaled = maxvalDFF .*(Phi * x_values);
        dFFRescaled = maxvalDFF .* dFF{ii};
    else
        dFFRescaled = dFF{ii};
        dataReconstructedRescaled = (Phi * x_values);
    end
   
    
    flatDFF = reshape(dFFRescaled, numel(dFFRescaled),[]);
    flatReconstructed = reshape(dataReconstructedRescaled,numel(dataReconstructedRescaled),[]);
    rval = corrcoef(flatDFF, flatReconstructed);
    varExpl = (rval(1,2))^2;
    disp(varExpl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%