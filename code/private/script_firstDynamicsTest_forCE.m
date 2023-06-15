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
%% Set up pointers to the fluorescence data and meta behavioral data

dat.dir   = 'C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegans/Data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

% dat.dFF = getNPYarray([dat.dir, dat.Ffile]);                               % Time x Neuron matrix of fluorescent activity
% dat.bhv = getNPYarray([dat.dir, dat.Bfile]);                               % Time x 1 vector of behavioral data
% dat.vel = getNPYarray([dat.dir, dat.Vfile]);                               % Time x 1 vector of velocity

dat = load('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegansData/WT_Stim.mat');
dat2 = load('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegansData/WT_NoStim.mat');

dFF = dat.WT_Stim(1).traces;
bhv = dat.WT_Stim(1).States;

% dFF_S_2 = dat.WT_Stim(2).traces;
% bhv_S_2 = dat.WT_Stim(2).States;
% 
% dFF_S_3 = dat.WT_Stim(3).traces;
% bhv_S_3 = dat.WT_Stim(3).States;
% 
% dFF_S_4 = dat.WT_Stim(4).traces;
% bhv_S_4 = dat.WT_Stim(4).States;
% 
% dFF_S_5 = dat.WT_Stim(5).traces;
% bhv_S_5 = dat.WT_Stim(5).States;
% 
% dFF_S_6 = dat.WT_Stim(6).traces;
% bhv_S_6 = dat.WT_Stim(6).States;
% 
% dFF_S_7 = dat.WT_Stim(7).traces;
% bhv_S_7 = dat.WT_Stim(7).States;

dFF2 = dat2.WT_NoStim(1).traces;
bhv2 = dat2.WT_NoStim(1).States;

% dFF2_NS_2 = dat2.WT_NoStim(2).traces;
% bhv2_NS_2 = dat2.WT_NoStim(2).States;
% 
% dFF2_NS_3 = dat2.WT_NoStim(3).traces;
% bhv2_NS_3 = dat2.WT_NoStim(3).States;
% 
% dFF2_NS_4 = dat2.WT_NoStim(4).traces;
% bhv2_NS_4 = dat2.WT_NoStim(4).States;
% 
% dFF2_NS_5 = dat2.WT_NoStim(5).traces;
% bhv2_NS_5 = dat2.WT_NoStim(5).States;



% %%
% fig10 = figure(10);
% 
% perplexityvals = [10,50,100,200];
% datasets = {dFF, dFF2, dFF(:,[41 104]), dFF2(:,[36 37])};
% bhvsets = {bhv, bhv2, bhv, bhv2};
% trialtype = {'Stim all', 'NoStim all', 'Stim AVA', 'NoStim AVA'};
% 
% for i = 1:4
%     dataset = datasets{i};
%     bhvset = bhvsets{i};
%     
%     for j = 1:4
%         subplot(4,4,(i-1)*4+j);
%         Y_tsne = tsne(dataset, 'perplexity', perplexityvals(j));
%         gscatter(Y_tsne(:,1),Y_tsne(:,2), bhvset);
%         title(sprintf('%s, plxty %d', trialtype{i}, perplexityvals(j)));
%         box off;      
%         set(gca, 'TickDir', 'out');
%         lg.Box = 'off';
%         axis tight;
% 
%     end
%   
% end
% %%
% % Connor Meehan, Jonathan Ebrahimian, Wayne Moore, and Stephen Meehan 
% % (2021). Uniform Manifold Approximation and Projection (UMAP) 
% % (https://www.mathworks.com/matlabcentral/fileexchange/71902), 
% % MATLAB Central File Exchange.
% 
% dFF_for_umap = [dFF bhv.'];
% dFF2_for_umap = [dFF2 bhv2.'];
% dFF_for_umap_AVA = [dFF(:,[41 104]) bhv.'];
% dFF2_for_umap_AVA = [dFF2(:,[36 37]) bhv2.'];
% 
% % fig11 = figure(11);
% 
% spreadvals = [0.1, 0.2, 0.3, 0.5];
% datasets = {dFF_for_umap, dFF2_for_umap, dFF_for_umap_AVA, dFF2_for_umap_AVA};
% trialtype = {'Stim all', 'NoStim all', 'Stim AVA', 'NoStim AVA'};
% 
% for i = 1:4
%     dataset = datasets{i};
%     
%     for j = 1:4
% %         subplot(4,4,(i-1)*4+j);
%         [reduction_out, umap_out, clusterIdentifiers_out, extras_out]=run_umap(dataset, 'label_column', 'end', 'min_dist', 0.1, 'spread', spreadvals(j), 'plot_title', sprintf('%s, spread %.2f', trialtype{i}, spreadvals(j))); 
% %         title(sprintf('%s, spread %.2f', trialtype{i}, spreadvals(j)));
% %         box off;      
% %         set(gca, 'TickDir', 'out');
% %         lg.Box = 'off';
% %         axis tight;
% 
%     end
%   
% end





% fig7 = figure(7);
% histogram(dFF);
% saveas(fig7, "fig7_ce_hist.png");
% 
% fig8 = figure(8);
% histogram(real(log10(dFF)));
% saveas(fig8, "fig8_ce_hist_log10.png");

%%
% dFF(abs(dFF)<1e-4 | dFF<0)=0; % remove small and negative values

% disp('no stim case');
% dFF = dFF2;

dFF(dFF<0)=0; % remove negative values


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

figure(2)
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
% %% Basic dimensionality reduction
% 
% tmpCentered = bsxfun(@plus, dFF, -mean(dFF,2));
% tmpCentered = tmpCentered.'*tmpCentered;  
% [U, D] = eigs(double(tmpCentered), 100); % reduced  from 200 to 100 for C Elegans
% 
% %%
% plot(0:(numel(diag(D))),[0;cumsum(diag(D))./(trace(tmpCentered))])
% set(gca, 'XLim',[0,numel(diag(D))+1])
% xlabel('PCA components included')
% ylabel('Variance Explained')
% box off

% clear tmpCentered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Set parameters
inf_opts.nF              = 10;
inf_opts.N               = 10; %scaled down from 100 (zebrafish, whole brain)
inf_opts.M               = size(dFF,2) ;
inf_opts.lambda_val      = 0.1; % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = 0.90;
inf_opts.lambda_b        = 0.01;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.max_iter2       = 10; %500 %why is this max_iter2 instead of max_iters? Why is max_iters 3000 if I don't set it here?
inf_opts.max_iters       = 100;
inf_opts.special = '';
% inf_opts.special = 'noobs';


%%
% partition - 5 fold cross-validation, 1st Stim Worm

dFF(dFF<0)=0; % remove negative values

rng('default')
c = cvpartition(size(dFF,1),'KFold',5);

for i = 1:5
    trainset = dFF(training(c,i),:);
    testset = dFF(test(c,i),:);
    
    [Phi, F] = bpdndf_dynamics_learning(trainset.', [], [], inf_opts);              % Run the dictionary learning algorithm

    sizedFF = size(testset.');
    dFFcell = mat2cell(testset.', sizedFF(1), sizedFF(2));
    [~,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients

%     c_values = B_cell{1}.';
%     x_values = A_cell{1};
%     x_values_timebyx = x_values.';
%     dataReconstructed = (Phi * x_values).';

    x_vec = create_reco(F, B_cell{1}, dFFcell{1}(:,1), 'cum', dFFcell{1},0);
    dataReconstructed = x_vec.';
    
    flatDFF = reshape(testset, numel(testset),[]);
    flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
    resid = flatDFF - flatReconstructed;
    SSE = sum(resid.^2);
    diffYFromMean = flatDFF - mean(flatDFF);
    SST = sum(diffYFromMean.^2);
    varExpl = 1-(SSE/SST);
    disp(varExpl);

%     %% post-processing - credit to Noga Mudrik for post_proc_main.m and visualize_dyns.m
%     c_vals_cell = {c_values.'};
%     ce_struct = struct('D',Phi.','F_cell', {F},'c_i',{c_vals_cell},'x', x_values, 'true_mat', testset.', 'app_mat', dataReconstructed.');
%     [vargout, app_mat, true_mat, additional_return]  = post_proc_main(ce_struct, [], 'post_processing_results_ce_3000iters_smalllambda.mat');
    save(sprintf('saveCEWorkspace_%diters_fold%d.mat',inf_opts.max_iters,i))
end
    


%%
% do dynamics transfer from worm to worm? (Stim case)

combos = nchoosek(1:7, 2);
for comboIndex = 1:7
    firstWormIndex = combos(comboIndex, 1);
    secondWormIndex = combos(comboIndex, 2); 
    
    trainset = dat.WT_Stim(firstWormIndex).traces;
    testset = dat.WT_Stim(secondWormIndex).traces;
    
    trainset(trainset<0)=0; % remove negative values
    testset(testset<0)=0; % remove negative values
    
    [Phi, F] = bpdndf_dynamics_learning(trainset.', [], [], inf_opts);              % Run the dictionary learning algorithm

    sizedFF = size(testset.');
    dFFcell = mat2cell(testset.', sizedFF(1), sizedFF(2));
    [~,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients

%     c_values = B_cell{1}.';
%     x_values = A_cell{1};
%     x_values_timebyx = x_values.';
%     dataReconstructed = (Phi * x_values).';
%     
    x_vec = create_reco(F, B_cell{1}, dFFcell{1}(:,1), 'cum', dFFcell{1},0);
    dataReconstructed = x_vec.';
    
    flatDFF = reshape(testset, numel(testset),[]);
    flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
    resid = flatDFF - flatReconstructed;
    SSE = sum(resid.^2);
    diffYFromMean = flatDFF - mean(flatDFF);
    SST = sum(diffYFromMean.^2);
    varExpl = 1-(SSE/SST);
    disp(varExpl);

%     %% post-processing - credit to Noga Mudrik for post_proc_main.m and visualize_dyns.m
%     c_vals_cell = {c_values.'};
%     ce_struct = struct('D',Phi.','F_cell', {F},'c_i',{c_vals_cell},'x', x_values, 'true_mat', testset.', 'app_mat', dataReconstructed.');
%     [vargout, app_mat, true_mat, additional_return]  = post_proc_main(ce_struct, [], 'post_processing_results_ce_3000iters_smalllambda.mat');
    save(sprintf('saveCEWorkspace_Stim_%diters_trainW%d_testW%d.mat',inf_opts.max_iters,firstWormIndex,secondWormIndex))
end
    
%%
% [Phi, F, allPhi, allF] = bpdndf_dynamics_learning(dFF.', [], [], inf_opts);              % Run the dictionary learning algorithm

[Phi, F] = bpdndf_dynamics_learning(dFF.', [], [], inf_opts);              % Run the dictionary learning algorithm

%%
sizedFF = size(dFF.');
dFFcell = mat2cell(dFF.', sizedFF(1), sizedFF(2));
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%%
figure(4)
cvalues_timebyc = B_cell{1}.';
plot(cvalues_timebyc) % plot c values over time
legend
%%
fig11 = figure(11);

perplexityvals = [10,50,100,200];
datasets = {cvalues_timebyc};
bhvsets = {bhv};
trialtype = {'Stim all'};

for i = 1
    dataset = datasets{i};
    bhvset = bhvsets{i};
    
    for j = 1:4
        subplot(1,4,(i-1)*4+j);
        Y_tsne = tsne(dataset, 'perplexity', perplexityvals(j));
        gscatter(Y_tsne(:,1),Y_tsne(:,2), bhvset);
        title(sprintf('%s, plxty %d', trialtype{i}, perplexityvals(j)));
        box off;      
        set(gca, 'TickDir', 'out');
        lg.Box = 'off';
        axis tight;

    end
  
end
%%

%% FIXME: how to average across cross-validated dynamics models


% COMPARE TO LINDERMAN - activ. vs. neurons
nToPlot  = 10; % 20 for zebrafish
ixToPlot = randsample(size(Phi,2), nToPlot);
figure(5)
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

%%
c_values = B_cell{1}.';
x_values = A_cell{1};
x_values_timebyx = x_values.';
dataReconstructed = (Phi * x_values).';

%FIXME
flatDFF = reshape(dFF, numel(dFF),[]);
flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
resid = flatDFF - flatReconstructed;
SSE = sum(resid.^2);
diffYFromMean = flatDFF - mean(flatDFF);
SST = sum(diffYFromMean.^2);
varExpl = 1-(SSE/SST);
disp(varExpl);
%%
% video
nFrames = size(dFF,1);
figVideo = figure(6);

vid = VideoWriter('Celegans_WT_Stim_01_noneg_keepsmall.avi');
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
    bar(dFF(frame_i,:));
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
    scatter(frame_i-19:frame_i, bhv(frame_i-19:frame_i));
    xticklabels({'-20','-15','-10','-5','t'});
    xlabel('time');
    axis([frame_i-19 frame_i 0 4]);
    yticks([0 1 2 3 4]);
    yticklabels({'0','1','2','3','4'});
    ylabel('State (Kato 2015)');
    title('Behavior');
    set(gca, 'TickDir', 'out');
    box off;
    
    writeVideo(vid, getframe(figVideo));
    
%     drawnow
end

close(vid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%