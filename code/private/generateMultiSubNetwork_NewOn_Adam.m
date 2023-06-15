%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data
% use 10 latent states, 10 "recorded channels", 3000 timepoints

dimSS = 3;
simulatedDmatrix = eye(2*dimSS,2*dimSS);

[latentStatesX, Fgt] = generateMultiSubNetwork_NewOn(1500,3,3,10,false,false,true);
dFF = cell(numel(latentStatesX),1);
for ll = 1:numel(latentStatesX)
    dFF{ll} = (simulatedDmatrix*latentStatesX{ll}).' + 0.0001*randn(size(latentStatesX{ll},2), size(simulatedDmatrix,1));
    dFF{ll} = dFF{ll}.';
end

figure(3); 
for ll = 1:5
    subplot(5,1,ll), plot([sum(dFF{ll}(1:dimSS,:).^2,1);sum(dFF{ll}(dimSS+(1:dimSS),:).^2,1) ].')
end
%%
%dFF(dFF<0)=0; % remove negative values - this is for fluorescence data
% dFF = dFF./max(abs(dFF),[],'all'); %set to abs
Fgt2 = cell(numel(Fgt)+1,1);
for ll = 1:3; Fgt2{ll} = zeros(2*dimSS,2*dimSS); end
Fgt2{1}(1:dimSS,1:dimSS)                 = Fgt{1};
Fgt2{2}(dimSS+(1:dimSS),dimSS+(1:dimSS)) = Fgt{2};
Fgt2{3}(dimSS+(1:dimSS),dimSS+(1:dimSS)) = Fgt{3};
Fgt2{4}            = 0.001*randn(2*dimSS,2*dimSS);

%%
figure(), cla;                            
subplot(2,3,[1,2,3]), imagesc(dFF{randsample(numel(dFF),1)})  
title('Simulated data dFF (scaled by max abs)')
colorbar;
xlabel('time points')
ylabel('channels')
subplot(2,3,4), imagesc(Fgt{1})
title('f1stdzd - top left of F')
subplot(2,3,5), imagesc(Fgt{2})
title('f2stdzd - bottom right, time 1 to T/2')
subplot(2,3,6), imagesc(Fgt{3})
title('f3stdzd - bottom right, time T/2+1 to T')
%%

% Set parameters
inf_opts.nF              = 4; % 4
inf_opts.N               = 2*dimSS; %scaled down from 100 (zebrafish, whole brain)
inf_opts.M               = size(dFF,2) ;
% disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = 0.0001;        % ASC added 7/25 % VARY - tradeoff: approx. SNR - see ratios of traces % 0.0001 (don't shrink this further - gets to be essentially 0 for solver_L1RLS)
inf_opts.lambda_history  = 2.2;           % ASC added 7/25
inf_opts.lambda_b        = 0.9;           % ASC added 7/25 %0.025
inf_opts.lambda_historyb = 1.5;           % ASC added 7/25 %0.7
inf_opts.tol             = 1e-8;          % 1e-3
inf_opts.max_iter2       = 20;            %500 %20
inf_opts.max_iters       = 3000;
inf_opts.special         = '';
inf_opts.F_update        = true;
inf_opts.D_update        = false;
inf_opts.N_ex            = 30;    % 40 %100 %30
inf_opts.T_s             = 200; % 30 %100 %20

inf_opts.step_d          = 1;  % ASC added 7/25
inf_opts.step_f          = 0.1;  % ASC added 7/25 % 30
inf_opts.plot_option     = 10;  % ASC added 7/25
inf_opts.lambda_f        = 10e-1;

inf_opts.solver_type    = 'tfocs'; %default: 'tfocs'
inf_opts.CVX_Precision  = 'low';   %default = 'default'
inf_opts.special        = 'noobs';

%%

[Phi, F] = bpdndf_dynamics_learning(dFF, [], eye(2*dimSS), inf_opts);              % Run the dictionary learning algorithm
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%
tSel = 2;
figure(1), cla; 
subplot(4,4,[1,4]), plot(dFF{tSel}.')   % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('Signal (GT)')
subplot(4,4,[5,8]), plot(A_cell{tSel}.')   % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('Reconstruction')
subplot(4,4,[9,12]), plot(B_cell{tSel}.')   % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('Dynamics coefficient amplitude')
legend('C1','C2','C3','C4')
subplot(4,4,13), imagesc(F{1}); axis image; axis off; colormap gray
subplot(4,4,14), imagesc(F{2}); axis image; axis off; colormap gray
subplot(4,4,15), imagesc(F{3}); axis image; axis off; colormap gray
subplot(4,4,16), imagesc(F{4}); axis image; axis off; colormap gray

%%
[A_cell2,B_cell2] = parallel_bilinear_dynamic_inference({(dFF.')*10}, simulatedDmatrix, Fgt2(1:inf_opts.nF), ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%
figure(2), cla;
subplot(3,4,[1,4]), imagesc(A_cell2{1})   % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('Signal Amplitude')
subplot(3,4,[5,8]), plot(B_cell2{1}.')   % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('Dynamics coefficient amplitude')
legend('C1','C2','C3','C4')
subplot(3,4,9), imagesc(Fgt2{1}); axis image; axis off; colormap gray
subplot(3,4,10), imagesc(Fgt2{2}); axis image; axis off; colormap gray
subplot(3,4,11), imagesc(Fgt2{3}); axis image; axis off; colormap gray
subplot(3,4,12), imagesc(Fgt2{4}); axis image; axis off; colormap gray

%%

nTimePoints = size(latentStatesX,2);
figure(), cla;
subplot(3,1,1);
xNormLatent = zeros(nTimePoints,1);
for i = 1:nTimePoints
    xNormLatent(i) = norm(latentStatesX(:,i));
end
plot(xNormLatent);
title('Norm(latent X)');
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
title('Norm(inferred X)');
xlabel('time points')
% colorbar;

%%

kTest = 20;
fDyn = [];
for ll = 1:numel(Fgt2)
    fDyn(:,ll) = Fgt2{ll}*(dFF(kTest,:).');
end

FISTAopts.pos        = false;
FISTAopts.lambda     = 0.1;
FISTAopts.check_grad = 0;
res = fista_lasso((dFF(kTest+1,:).'), fDyn, [], FISTAopts)
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
