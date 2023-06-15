%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code

% addpath(genpath('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/MIND_example_for_adam/for_adam/helpers/'));
addpath(genpath('.'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

% dat.dFF = getNPYarray([dat.dir, dat.Ffile]);                               % Time x Neuron matrix of fluorescent activity
% dat.bhv = getNPYarray([dat.dir, dat.Bfile]);                               % Time x 1 vector of behavioral data
% dat.vel = getNPYarray([dat.dir, dat.Vfile]);                               % Time x 1 vector of velocity

% from embed_lorenz (Low et al. 2018)
[X Y Z] = lorenz(28, 10, 8/3, [0 1 1.05]);

% dI = 10;
% double the sampling rate?
dI = 5;
dFF = zeros(length(X(1:dI:end)),3);
dFF(:,1) = X(1:dI:end);
dFF(:,2) = Y(1:dI:end);
dFF(:,3) = Z(1:dI:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot some data
fig1 = figure(1);
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

saveas(fig1, "fig1_lorenz_upsampled_datamatrix.png");


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fig2 = figure(2);
% subplot(1,2,1); plot(bhv)
% box off
% title('Behavioral data')
% set(gca,'XLim', [1, numel(bhv)])
% 
%  subplot(1,2,2); 

plot(dFF);
% subplot(1,2,2); plot(bsxfun(@plus, 0:2:6, 10*dat.dFF(randsample(1:4,4),:)'))
box off
xlabel('Frame number')
ylabel('DF/F')
title('Example Traces')
% set(gca,'XLim', [1, size(dat.dFF,2)])
set(gca,'XLim', [1, size(dFF,1)])

saveas(fig2, "fig2_lorenz_upsampled_dFF_exampletraces.png");


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
inf_opts.lambda_val      = 0.1; % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = 0.90;
inf_opts.lambda_b        = 0.01;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3;
inf_opts.max_iter2       = 10; %500 %why is this max_iter2 instead of max_iters? Why is max_iters 3000 if I don't set it here?
inf_opts.special = '';

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

saveas(fig4, "fig5_lorenz_upsampled_cvalues.png");

fig11 = figure(11);

perplexityvals = [10,50,100,200];
datasets = {cvalues_timebyc};

for i = 1
    dataset = datasets{i};
    
    for j = 1:4
        subplot(1,4,(i-1)*4+j);
        Y_tsne = tsne(dataset, 'perplexity', perplexityvals(j));
        gscatter(Y_tsne(:,1),Y_tsne(:,2));
        title(sprintf('plxty %d', perplexityvals(j)));
        box off;      
        set(gca, 'TickDir', 'out');
        lg.Box = 'off';
        axis tight;

    end
  
end

saveas(fig11, "fig11_lorenz_upsampled_tsne.png");
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

saveas(fig5, "fig5_lorenz_upsampled.png");


%%
c_values = B_cell{1}.';
x_values = A_cell{1};
x_values_timebyx = x_values.';
dataReconstructed = (Phi * x_values).';

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

save('saveLorenzWorkspace_upsampled.mat')

%%
% video


nFrames = size(dFF,1);

%nFrames = 500;
figVideo = figure(6);

vid = VideoWriter('Lorenz_upsampled.avi');
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
    
    % from embed_lorenz
    subplot(3,2,5);
    err = sqrt(sum((dataReconstructed-dFF).^2,2));
    scatter3(dFF(:,1), dFF(:,2), dFF(:,3),[],err,'.')
    title('Lorenz colored by reconstruction error')
    set(gca, 'TickDir', 'out');
    box off;
    
    % from embed_lorenz
    subplot(3,2,6);
    scatter3(dFF(frame_i-19:frame_i,1), dFF(frame_i-19:frame_i,2), dFF(frame_i-19:frame_i,3),[],err(frame_i-19:frame_i),'.')
    title('Lorenz: reconstruction error, last 20 frames')
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
imshow(normalizedDFF(1:100,:).');
title('Real Lorenz vals');

subplot(2,1,2)
normalizedRecon = zeros(size(dataReconstructed,1),size(dataReconstructed,2));
for i = 1:3
    normalizedRecon(:,i) = (dataReconstructed(:,i) - min(dataReconstructed)) / ( max(dataReconstructed) - min(dataReconstructed) );
end
imshow(normalizedRecon(1:100,:).');
title('Reconstructed');

saveas(fig8, "fig8_lorenz_upsampled_first100vals_realvsreconstructed.png");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%