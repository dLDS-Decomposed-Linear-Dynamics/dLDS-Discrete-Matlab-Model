%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('../../../../../../home/eyezere1/npy-matlab/npy-matlab/'))                % Use steinmetz code for loady npy into matlab
addpath(genpath('.'))
% addpath(genpath('../../../../home/adamsc/GITrepos/zebrafishDynamics/code/'))            % This is the main codebase for the project
% addpath(genpath('../../../../home/adamsc/GITrepos/dynamics_learning/'))                 % This is the package for learning dynamics 

addpath(genpath('../../../../../../home/adamsc/data/zebrafish'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data

dat.dir = '../../../../../../home/adamsc/data/zebrafish/';
dat.Ffile = 'for_JHU_cells_dff.npy';
dat.Bfile = 'for_JHU_gad1b-6hz-backwardpulse_integrate_behavior_ds.npy';
dat.Vfile = 'for_JHU_-gad1b-6hz-backwardpulse_integrate_visual_velocity_ds.npy';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

dat.dFF = readNPY([dat.dir, dat.Ffile]);   % getNPYarray doesn't exist % Time x Neuron matrix of fluorescent activity
dat.bhv = readNPY([dat.dir, dat.Bfile]);                               % Time x 1 vector of behavioral data
dat.vel = readNPY([dat.dir, dat.Vfile]);                               % Time x 1 vector of velocity

fig7 = figure(7);
histogram(dat.dFF);
saveas(fig7, "fig7_zebrafish_hist.png");

fig8 = figure(8);
histogram(real(log10(dat.dFF)));
saveas(fig8, "fig8_zebrafish_hist_log10.png");

fprintf('number of zero values %f\n', sum(dat.dFF(:)==0));
fprintf('percentage of zero values %f\n', sum(dat.dFF(:)==0)/numel(dat.dFF));
% minus3sigma = mean(real(log10(dat.dFF(:))))-std(real(log10(dat.dFF(:))));
% fprintf('number of values to set to 0 %f\n', sum(abs(dat.dFF(:))< minus3sigma | dat.dFF(:)<0));
fprintf('number of values to set to 0 %f\n', sum(dat.dFF(:)<0));
% dat.dFF(abs(dat.dFF(:))< minus3sigma | dat.dFF(:)<0)=0; % remove small and negative values
dat.dFF(dat.dFF(:)<0)=0; % remove negative values
fprintf('total number of values %f\n', numel(dat.dFF));
fprintf('percentage of zero values %f\n', sum(dat.dFF(:)==0)/numel(dat.dFF));

fig9 = figure(9);
histogram(real(log10(dat.dFF)));
saveas(fig9, "fig9_zebrafish_hist_log10_noneg_keepsmall.png");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot some data
fig1 = figure(1);
subplot(2,1,1), imagesc(dat.dFF, [0,1])
axis image; box off
title('Data Matrix')
xlabel('Time? (frames)')
ylabel('Neuron ID?')


subplot(2,1,2), imagesc(diag(1./max(dat.dFF,[],2))*dat.dFF, [0,2])
axis image; box off
title('Normalized Data Matrix')
xlabel('Time? (frames)')
ylabel('Neuron ID?')

saveas(fig1, "fig1_zebrafish_noneg_keepsmall.png");

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fig2 = figure(2);
subplot(1,2,1); plot(dat.bhv)
box off
title('Behavioral data')
set(gca,'XLim', [1, numel(dat.bhv)])

subplot(1,2,2); plot(bsxfun(@plus, 0:2:6, 10*dat.dFF(randsample(1:4,4),:)'))
box off
xlabel('Frame number')
ylabel('DF/F')
title('Example Traces')
set(gca,'XLim', [1, size(dat.dFF,2)])

saveas(fig2, "fig2_zebrafish_noneg_keepsmall.png");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic dimensionality reduction

tmpCentered = bsxfun(@plus, dat.dFF, -mean(dat.dFF,2));
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
inf_opts.N               = 100;
inf_opts.M               = size(dat.dFF,2) ;
inf_opts.lambda_val      = 0.1;
inf_opts.lambda_history  = 0.90;
inf_opts.lambda_b        = 0.01;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.max_iter2       = 500;
inf_opts.special = '';

[Phi, F] = bpdndf_dynamics_learning(dat.dFF.', [], [], inf_opts);              % Run the dictionary learning algorithm

% %%
% 
% nToPlot  = 20;
% ixToPlot = randsample(size(Phi,2), nToPlot);
% fig4 = figure(4);
% subplot(2,1,1); plot(bsxfun(@plus, -Phi(:,ixToPlot)*diag(1./max(Phi(:,ixToPlot),[],1)), (1:nToPlot)-1))
% box off
% set(gca, 'XLim', [1,size(Phi,1)])
% xlabel('Neuron index')
% ylabel('Rel. contribution')
% title('Example activation profiles')
% 
% Ftmp = [];
% for ll = 1:numel(F); Ftmp = cat(2,Ftmp, reshape(F{ll},[],1)); end
% 
% subplot(2,1,2), imagesc(basis2img2(Ftmp, size(F{1}), [1,5]))
% axis image
% axis off
% title('Learned dynamics functions')
% clear nToPlot ixToPlot Ftmp
% saveas(fig4, "fig4_zebrafish.png");
%%
sizedFF = size(dat.dFF.');
dFFcell = mat2cell(dat.dFF.', sizedFF(1), sizedFF(2));
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%%
fig4 = figure(4);
plot(B_cell{1}.') % plot c values over time
legend
saveas(fig4, "fig4_zebrafish_noneg_keepsmall.png");

% COMPARE TO LINDERMAN - activ. vs. neurons
nToPlot  = 10; % 20 for zebrafish
ixToPlot = randsample(size(Phi,2), nToPlot);
fig5 = figure(5);
subplot(2,1,1); 
% plot(bsxfun(@plus, -Phi(:,ixToPlot)*diag(1./max(Phi(:,ixToPlot),[],1)), (1:nToPlot)-1))
plot(bsxfun(@plus, ...
    Phi(:,ixToPlot)*diag(1./max(abs(Phi(:,ixToPlot)),[],1)), (1:nToPlot)-1))


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
saveas(fig5, "fig5_zebrafish_noneg_keepsmall.png");

%%
% video
c_values = B_cell{1}.';
x_values = A_cell{1};
x_values_timebyx = x_values.';
dataReconstructed = (Phi * x_values).';


save('saveZebrafishWorkspace.mat')

nFrames = size(dat.dFF,1);
figVideo = figure(6);

vid = VideoWriter('Zebra_for_JHU_cells_dff_noneg_keepsmall.avi');
vid.Quality = 100;
vid.FrameRate = 10;
open(vid);
for frame_i = 20:nFrames
    subplot(3,2,1);
    bar(dataReconstructed(frame_i,:));
    xlabel('Neuron ID');
    ylabel('Activation');
    title('Reconstructed');
    set(gca, 'TickDir', 'out');
    box off;
    
    subplot(3,2,3);
    bar(dat.dFF(frame_i,:));
    xlabel('Neuron ID');
    ylabel('Activation');
    title('Recorded');
    set(gca, 'TickDir', 'out');
    box off;
    
    subplot(3,2,2);
    plot(c_values(frame_i-19:frame_i,:)); 
    xticklabels({'-20','-15','-10','-5','t'});
    xlabel('time');
    title('c values');
    set(gca, 'TickDir', 'out');
    box off;
    
    subplot(3,2,4);
    plot(x_values_timebyx(frame_i-19:frame_i,:)); 
    xticklabels({'-20','-15','-10','-5','t'});
    xlabel('time');
    title('x values');
    set(gca, 'TickDir', 'out');
    box off;
    
    subplot(3,2,5);
    % zebrafishBhv = dat.bhv(1,frame_i-19:frame_i);
    zebrafishBhv = dat.bhv(frame_i-19:frame_i);
    plot(zebrafishBhv/std(zebrafishBhv))
    xticklabels({'-20','-15','-10','-5','t'});
    xlabel('time');
    title('Behavior trace (bhv/std(bhv))');
    set(gca, 'TickDir', 'out');
    box off;
    
    % subplot(3,2,6);
    % zebrafishDrift = dat.bhv(2,frame_i-19:frame_i);
    % plot(zebrafishDrift); 
    % xticklabels({'-20','-15','-10','-5','t'});
    % xlabel('time');
    % title('Drift pulses');
    % set(gca, 'TickDir', 'out');
    % box off;
    
    writeVideo(vid, getframe(figVideo));
    
%     drawnow
end

close(vid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%