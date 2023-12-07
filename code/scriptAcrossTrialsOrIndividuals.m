%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data

% dat.dir   = 'C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegans/Data';

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
% inf_opts.AcrossIndividuals     = false;

% 1 worm, multiple durations of trials (artificially trimmed) (same channels, could be different trial durations)
% whichStimWorm = 1;
% dFF = {dFF_S_AllWorms{whichStimWorm}(:,1:2000);dFF_S_AllWorms{whichStimWorm}(:,1:1000);dFF_S_AllWorms{whichStimWorm}(:,1:500)}; 
% bhv = {bhv_S_AllWorms{whichStimWorm}(:,1:2000);bhv_S_AllWorms{whichStimWorm}(:,1:1000);bhv_S_AllWorms{whichStimWorm}(:,1:500)}; 
% inf_opts.AcrossIndividuals     = false;

% % multiple worms (different channels, different trial durations)
dFF = dFF_S_AllWorms;
bhv = bhv_S_AllWorms;
inf_opts.AcrossIndividuals     = true;

%%
for ii = 1:size(dFF,1)
    thisDFF = dFF{ii};
    thisDFF(thisDFF<0) = 0;
    thisDFF = thisDFF./max(thisDFF,[],'all');
    dFF{ii} = thisDFF;
end


%%

% Set parameters
inf_opts.nF              = input('nF:'); % 10 per worm
inf_opts.N               = input('n:'); % 10 
disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = input('lambda_val:'); % 0.1 % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = input('lambda_history:'); % 0.9
inf_opts.lambda_b        = input('lambda_b:'); % 0.01 (also other reg params) - 0.1 all disappear, 0.05 only 2 traces
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3; % 1e-3
inf_opts.max_iter2       = 50; %50 
inf_opts.max_iters       = input('max_iters:'); %500
inf_opts.special         = '';
inf_opts.F_update        = true;
inf_opts.D_update        = true;
% inf_opts.AcrossIndividuals     = false;
inf_opts.N_ex            = input('N_ex:');%40;
inf_opts.T_s             = input('T_s:');%30;
% inf_opts.sampleProportionally = true; %default: prop by trial length
% inf_opts.special         = 'noobs';

% inf_opts.step_d          = 20;  % ASC added 7/25
% inf_opts.step_f          = 30;  % ASC added 7/25
% inf_opts.plot_option     = 20;  % ASC added 7/25
% inf_opts.lambda_f        = 0.05; % ASC added 7/28


if inf_opts.AcrossIndividuals
    [inf_opts.M,~]               = cellfun(@size,dFF,'UniformOutput',false);%size(dFF,1) ; % EY changed 9/6/23 for ManyWorms case
else
    inf_opts.M = size(dFF{1},1) ; % original dimension (# channels)
end

%%

[Phi, F] = bpdndf_dynamics_learning(dFF, [], [], inf_opts);              % Run the dictionary learning algorithm


if inf_opts.AcrossIndividuals
    for ii = 1:size(Phi,1)
        [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference(dFF(ii), Phi{ii}, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
    end
else
    [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients

end

%%
if inf_opts.AcrossIndividuals
    for ii = 1:size(Phi,1)
        x_values = cell2mat(A_cell{ii});
        dataReconstructed = (Phi{ii} * x_values).';    
        
        flatDFF = reshape(dFF{ii}, numel(dFF{ii}),[]);
        flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
        rval = corrcoef(dFF{ii}, dataReconstructed.');
        varExpl = (rval(1,2))^2;
        disp(varExpl);

    end
else
    for ii = 1:size(A_cell,1) %in the event of multiple trials
        x_values = A_cell{ii};
        dataReconstructed = (Phi * x_values).';    
        
        maxvalDFF = max(dFF{ii},[],'all');
        stdzdOn = false;
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
        % pause();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%