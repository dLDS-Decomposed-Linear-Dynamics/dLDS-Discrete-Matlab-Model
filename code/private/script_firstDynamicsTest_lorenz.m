%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code

addpath(genpath('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/MIND_example_for_adam/for_adam/helpers/'));
addpath(genpath('.'))
addpath(genpath('../..'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

% dat.dFF = getNPYarray([dat.dir, dat.Ffile]);                               % Time x Neuron matrix of fluorescent activity
% dat.bhv = getNPYarray([dat.dir, dat.Bfile]);                               % Time x 1 vector of behavioral data
% dat.vel = getNPYarray([dat.dir, dat.Vfile]);                               % Time x 1 vector of velocity

% from embed_lorenz (Low et al. 2018)
[X Y Z] = lorenz(28, 10, 8/3, [0 1 1.05]);

dI = 10;
dFF = zeros(length(X(1:dI:end)),3);
dFF(:,1) = X(1:dI:end);
dFF(:,2) = Y(1:dI:end);
dFF(:,3) = Z(1:dI:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot some data
% fig1 = figure(1);
% subplot(2,1,1), imagesc(dFF, [0,1])
% axis image; box off
% title('Data Matrix')
% xlabel('Time? (frames)')
% ylabel('Neuron ID?')
% 
% 
% subplot(2,1,2), imagesc(diag(1./max(dFF,[],2))*dFF, [0,2])
% axis image; box off
% title('Normalized Data Matrix')
% xlabel('Time? (frames)')
% ylabel('Neuron ID?')
% 
% saveas(fig1, "fig1_lorenz_datamatrix_3000iters_smalllambda_yequalsx.png");
% 
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
% 
% fig2 = figure(2);
% % subplot(1,2,1); plot(bhv)
% % box off
% % title('Behavioral data')
% % set(gca,'XLim', [1, numel(bhv)])
% % 
% %  subplot(1,2,2); 
% 
% plot(dFF);
% % subplot(1,2,2); plot(bsxfun(@plus, 0:2:6, 10*dat.dFF(randsample(1:4,4),:)'))
% box off
% xlabel('Frame number')
% ylabel('DF/F')
% title('Example Traces')
% % set(gca,'XLim', [1, size(dat.dFF,2)])
% set(gca,'XLim', [1, size(dFF,1)])
% 
% saveas(fig2, "fig2_lorenz_dFF_exampletraces_3000iters_smalllambda_yequalsx.png");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Basic dimensionality reduction
% 
% tmpCentered = bsxfun(@plus, dFF, -mean(dFF,2));
% tmpCentered = tmpCentered.'*tmpCentered;  
% [U, D] = eigs(double(tmpCentered), 3); % reduced  from 200 to 100 for C Elegans
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
inf_opts.nF              = 3;
inf_opts.N               = 3; %scaled down from 100 (zebrafish, whole brain)
inf_opts.M               = size(dFF,2) ;
inf_opts.lambda_val      = 1e-5; %was 0.1 % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = 1e-5; %was 0.90
inf_opts.lambda_b        = 0.01; %was 0.01, then 0.1 
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.max_iter2       = 10; %500 %why is this max_iter2 instead of max_iters? Why is max_iters 3000 if I don't set it here?
% inf_opts.max_iters       = 3000;
inf_opts.special = '';
inf_opts.type_solve      = 'lasso';
inf_opts.D_update        = false; 
inf_opts.T_s             = size(dFF,1); % size(dFF,1)% changes anyway to size of data (number of timepoints)(no it doesn't)
inf_opts.fix_f           = false;
inf_opts.step_f          = 0.001; % originally 20, then 1, then 0.01 (still explodes), then 0.001 - I think it hit local minimum

dFF = bsxfun(@plus, dFF, -mean(dFF,1));

% why are the c values exploding? See infer_B
[Phi, F, allPhi, allF] = bpdndf_dynamics_learning(dFF.', [], [], inf_opts);              % Run the dictionary learning algorithm

%%
sizedFF = size(dFF.');
dFFcell = mat2cell(dFF.', sizedFF(1), sizedFF(2));
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%%
fig4 = figure(4);
cvalues_timebyc = B_cell{1}.';
plot(cvalues_timebyc) % plot c values over time
title('c values by time');
legend

saveas(fig4, "fig4_lorenz_cvalues_3000iters_smalllambda_yequalsx_lambdaval0_lambdahistory0.png");

% fig11 = figure(11);
% 
% perplexityvals = [10,50,100,200];
% datasets = {cvalues_timebyc};
% 
% for i = 1
%     dataset = datasets{i};
%     
%     for j = 1:4
%         subplot(1,4,(i-1)*4+j);
%         Y_tsne = tsne(dataset, 'perplexity', perplexityvals(j));
%         gscatter(Y_tsne(:,1),Y_tsne(:,2));
%         title(sprintf('plxty %d', perplexityvals(j)));
%         box off;      
%         set(gca, 'TickDir', 'out');
%         lg.Box = 'off';
%         axis tight;
% 
%     end
%   
% end
% 
% saveas(fig11, "fig11_lorenz_tsne.png");
%%
% COMPARE TO LINDERMAN - activ. vs. neurons
nToPlot  = 3; % 20 for zebrafish
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

saveas(fig5, "fig5_lorenz_3000iters_smalllambda_yequalsx_lambdaval0_lambdahistory0.png");


%%
c_values = B_cell{1}.';
x_values = A_cell{1};
if iscell(x_values)
    x_values = squeeze(cell2mat(x_values));
end
% x_values_timebyx = x_values.';
% dataReconstructed = (Phi * x_values).';

%% Do reconstruction using G and c 
x_vec = create_reco(F, B_cell{1}, dFFcell{1}(:,1), 'cum', dFFcell{1},0);

dataReconstructed = x_vec.';
%%
flatDFF = reshape(dFF, numel(dFF),[]);
flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
resid = flatDFF - flatReconstructed;
SSE = sum(resid.^2);
diffYFromMean = flatDFF - mean(flatDFF);
SST = sum(diffYFromMean.^2);
varExpl = 1-(SSE/SST);
disp(varExpl);

% numeratorCorrCoeff = cov(dataReconstructed, dFF);
% denominatorCorrCoeff = std(dataReconstructed) * std(dFF);
% rCorrCoeff3D = numeratorCorrCoeff/denominatorCorrCoeff;
% combinedMean = sum(rCorrCoeff3D)/3; % same n in each col originally
% % combinedVar = (sum(var(dataReconstructed), var(dFF),))/3; %https://www.emathzone.com/tutorials/basic-statistics/combined-variance.html

save('saveLorenzWorkspace_3000iters_smalllambda_yequalsx_lambdaval0_lambdahistory0.mat')

%%
% video


nFrames = size(dFF,1)-1;

%nFrames = 500;
figVideo = figure(6);

vid = VideoWriter('Lorenz_3000iters_smalllambda_yequalsx_lambdaval0_lambdahistory0.avi');
vid.Quality = 100;
vid.FrameRate = 10;
open(vid);
for frame_i = 20:nFrames
    subplot(4,2,1);
    bar(dataReconstructed(frame_i,:));
    xlabel('Neuron ID');
    ylabel('Activation');
    title('Reconstructed');
    set(gca, 'TickDir', 'out');
    box off;
    
    subplot(4,2,3);
    bar(dFF(frame_i,:));
    xlabel('Neuron ID');
    ylabel('Activation');
    title('Recorded');
    set(gca, 'TickDir', 'out');
    box off;
    
    subplot(4,2,2);
    plot(c_values(frame_i-19:frame_i,:)); 
    xticklabels({'-20','-15','-10','-5','t'});
    xlabel('time');
    title('c values');
    set(gca, 'TickDir', 'out');
    box off;
    
    subplot(4,2,4);
    plot(x_values_timebyx(frame_i-19:frame_i,:)); 
    xticklabels({'-20','-15','-10','-5','t'});
    xlabel('time');
    title('x values');
    set(gca, 'TickDir', 'out');
    box off;
    
    % from embed_lorenz
    subplot(4,2,5);
    err = sqrt(sum((dataReconstructed-dFF).^2,2)); % reconstruction error
    scatter3(dFF(:,1), dFF(:,2), dFF(:,3),[],err,'.')
    title('Lorenz colored by reconstruction error')
    set(gca, 'TickDir', 'out');
    box off;
    
    % from embed_lorenz
    subplot(4,2,6);
    scatter3(dFF(frame_i-19:frame_i,1), dFF(frame_i-19:frame_i,2), dFF(frame_i-19:frame_i,3),[],err(frame_i-19:frame_i),'.')
    title('Lorenz: reconstruction error, last 20 frames')
    set(gca, 'TickDir', 'out');
    box off;
    colorbar;
    
    subplot(4,2,8);
    scatter3(dataReconstructed(frame_i-19:frame_i,1), dataReconstructed(frame_i-19:frame_i,2), dataReconstructed(frame_i-19:frame_i,3),[],err(frame_i-19:frame_i),'.')
    title('Reconstruction, last 20 frames')
    set(gca, 'TickDir', 'out');
    box off;
    colorbar;
    
    writeVideo(vid, getframe(figVideo));
    
%     drawnow
end


close(vid);

%% plot approximated values vs real over time
fig8 = figure(8);
subplot(2,1,1)
normalizedDFF = zeros(size(dFF,1),size(dFF,2));
for i = 1:3
    normalizedDFF(:,i) = (dFF(:,i) - min(dFF)) / ( max(dFF) - min(dFF) );
end
imagesc(normalizedDFF(1:100,:).');
title('Real Lorenz vals');

subplot(2,1,2)
normalizedRecon = zeros(size(dataReconstructed,1),size(dataReconstructed,2));
for i = 1:3
    normalizedRecon(:,i) = (dataReconstructed(:,i) - min(dataReconstructed)) / ( max(dataReconstructed) - min(dataReconstructed) );
end
imagesc(normalizedRecon(1:100,:).');
title('Reconstructed');

saveas(fig8, "fig8_lorenz_first100vals_realvsreconstructed_3000iters_smalllambda_yequalsx.png");


%% plot lorenz reconstruction itself
fig9 = figure(9);
err = sqrt(sum((dataReconstructed-dFF).^2,2));
scatter3(dataReconstructed(:,1), dataReconstructed(:,2), dataReconstructed(:,3),[],err,'.')
title('Lorenz reconstruction (color: error)')
set(gca, 'TickDir', 'out');
box off;
saveas(fig9, "fig9_lorenz_reconstructed_3000iters_smalllambda_yequalsx.png");

%% plot lorenz reconstruction itself
fig10 = figure(10);
scatter3(dataReconstructed(1:end-1,1), dataReconstructed(1:end-1,2), dataReconstructed(1:end-1,3),[],c_values/max(c_values),'.')
title('Lorenz reconstruction (color: c vals)')
set(gca, 'TickDir', 'out');
box off;
saveas(fig10, "fig10_lorenz_reconstructed_c_3000iters_smalllambda_yequalsx.png");

%% plot lorenz reconstruction itself
fig11 = figure(11);
[~,maxCIndexAtEachTimePt] = max(abs(c_values.'));
scatter3(dataReconstructed(1:end-1,1), dataReconstructed(1:end-1,2), dataReconstructed(1:end-1,3),[],maxCIndexAtEachTimePt,'.')
title('Lorenz reconstruction (color: argmax abs(c))')
set(gca, 'TickDir', 'out');
% legend 
box off;
saveas(fig11, "fig11_lorenz_reconstructed_argmaxc_3000iters_smalllambda_yequalsx.png");


%% post-processing - credit to Noga Mudrik for post_proc_main.m and visualize_dyns.m
% c_vals_cell = num2cell(c_values.');
c_vals_cell = {c_values.'};
lorenz_struct = struct('D',Phi.','F_cell', {F},'c_i',{c_vals_cell},'x', x_values, 'true_mat', dFF.', 'app_mat', dataReconstructed.');
[vargout, app_mat, true_mat, additional_return]  = post_proc_main(lorenz_struct, [], 'post_processing_results_lorenz_3000iters_smalllambda.mat');

scatter3(app_mat(:,1), app_mat(:,2), app_mat(:,3),[],[],'.')
title('Lorenz reconstruction')
set(gca, 'TickDir', 'out');
box off; 

% %% new reconstruction
% x_vec = create_reco(F, c_values, dFF(1,:), 'cum', [] ,0); %not sure what the ground truth dynamics would be here

%% plot G_k 
xVal = x_values(:,1);
applyDyns = zeros(size(F,2),10000,size(xVal,1));

%FIXME: set c_k as 1/||G_k|| to stabilize G_k
%FIXME: or keep the original cs
for dyn = 1:size(F,2)
    for i = 1:10000
        newXval = F{dyn}*xVal;
        applyDyns(dyn,i,:) = newXval;
        xVal = newXval;
    end  
    figure(20+dyn);
    scatter3(applyDyns(dyn,:,1), applyDyns(dyn,:,2), applyDyns(dyn,:,3));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%