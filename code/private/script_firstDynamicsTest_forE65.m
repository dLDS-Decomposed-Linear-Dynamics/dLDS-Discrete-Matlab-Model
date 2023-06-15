%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code

% addpath(genpath('/home/adam/GITrepos/npy-matlab/'))                        % Use steinmetz code for loady npy into matlab
% addpath(genpath('/home/adam/GITrepos/zebrafishDynamics/code/'))            % This is the main codebase for the project
% addpath(genpath('/home/adam/GITrepos/dynamics_learning/'))                 % This is the package for learning dynamics 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data

dat.dir   = 'C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/MIND_example_for_adam/for_adam';
% dat.Ffile = 'for_JHU_cells_dff.npy';
% dat.Bfile = 'for_JHU_gad1b-6hz-backwardpulse_integrate_behavior_ds.npy';
% dat.Vfile = 'for_JHU_-gad1b-6hz-backwardpulse_integrate_visual_velocity_ds.npy';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

% dat.dFF = getNPYarray([dat.dir, dat.Ffile]);                               % Time x Neuron matrix of fluorescent activity
% dat.bhv = getNPYarray([dat.dir, dat.Bfile]);                               % Time x 1 vector of behavioral data
% dat.vel = getNPYarray([dat.dir, dat.Vfile]);                               % Time x 1 vector of velocity

dat = load('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/MIND_example_for_adam/for_adam/E65.mat');

dFF = dat.nic_output.ROIactivities;
bhv = dat.nic_output.behavioralVariables;

assert(ismatrix(dFF));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot some data
figure(1)
subplot(2,1,1), imagesc(dFF, [0,1])
axis image; box off
title('Data Matrix')
xlabel('Time? (frames)')
ylabel('Neuron ID?')


subplot(2,1,2), imagesc(diag(1./max(dFF,[],2))*dFF, [0,2])
axis image; box off
title('Normalized Data Matrix')
xlabel('Time? (frames)')
ylabel('Neuron ID?')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% figure(2)
% subplot(1,2,1); plot(bhv)
% box off
% title('Behavioral data')
% set(gca,'XLim', [1, numel(bhv)])

 subplot(1,2,2); plot(bsxfun(@plus, 0:2:6, dFF(:,randsample(size(dFF,2),4))))
% subplot(1,2,2); plot(bsxfun(@plus, 0:2:6, 10*dat.dFF(randsample(1:4,4),:)'))
box off
xlabel('Frame number')
ylabel('DF/F')
title('Example Traces')
% set(gca,'XLim', [1, size(dat.dFF,2)])
set(gca,'XLim', [1, size(dFF,1)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic dimensionality reduction

tmpCentered = bsxfun(@plus, dFF, -mean(dFF,2));
tmpCentered = tmpCentered.'*tmpCentered;  
[U, D] = eigs(double(tmpCentered), 200);

%%
plot(0:(numel(diag(D))),[0;cumsum(diag(D))./(trace(tmpCentered))])
set(gca, 'XLim',[0,numel(diag(D))+1])
xlabel('PCA components included')
ylabel('Variance Explained')
box off

% clear tmpCentered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Set parameters
inf_opts.nF              = 5;
inf_opts.N               = 40; % scaled up bc approx. 4x as many neurons as C Elegans
inf_opts.M               = size(dFF,2) ;
inf_opts.lambda_val      = 0.1;
inf_opts.lambda_history  = 0.90;
inf_opts.lambda_b        = 0.01;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.max_iter2       = 500;
inf_opts.special = '';

[Phi, F] = bpdndf_dynamics_learning(dFF.', [], [], inf_opts);              % Run the dictionary learning algorithm

%%
sizedFF = size(dFF.');
dFFcell = mat2cell(dFF.', sizedFF(1), sizedFF(2));
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%%
figure(5)
plot(B_cell{1}.')
legend
%%

nToPlot  = 20;
ixToPlot = randsample(size(Phi,2), nToPlot);
figure(4)
subplot(2,1,1); plot(bsxfun(@plus, -Phi(:,ixToPlot)*diag(1./max(Phi(:,ixToPlot),[],1)), (1:nToPlot)-1))
box off
set(gca, 'XLim', [1,size(Phi,1)])
xlabel('Neuron index')
ylabel('Rel. contribution')
title('Example activation profiles')

Ftmp = [];
for ll = 1:numel(F); Ftmp = cat(2,Ftmp, reshape(F{ll},[],1)); end

subplot(2,1,2), imagesc(basis2img2(Ftmp, size(F{1}), [1,5]))
axis image
axis off
title('Learned dynamics functions')
clear nToPlot ixToPlot Ftmp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%