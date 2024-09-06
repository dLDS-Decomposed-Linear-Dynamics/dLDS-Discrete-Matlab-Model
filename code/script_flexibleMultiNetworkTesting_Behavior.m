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
%%
nBhv = 10;

% fourierPsi = fourierBasis(sum(nF,"all"));
% rowPsi = fourierPsi(rand);
% rowPsi = rowPsi(1,sum(nF,"all")+1:end);
% simulatedPsi = repmat(rowPsi,nBhv,1);

fourierPsi = fourierBasis(nBhv);
randomFrequencies = 10*rand(1,sum(nF,'all'));
simulatedPsi = zeros(nBhv,sum(nF,'all'));
for i = 1:sum(nF,'all')
    colPsi = fourierPsi(randomFrequencies(i)).';
    colPsi = colPsi(nBhv+1:end,1);
    simulatedPsi(:,i) = colPsi;
end

behaviorData = cellfun(@(x) simulatedPsi*x,groundTruthStates,'UniformOutput',false);
% disp('FIXME: can behaviorData be multiple cells too?')
dFF = statesToObservations(simulatedDmatrix, latentStatesX, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bayesopt

rng default
lambda_behavior = optimizableVariable('lambda_behavior',[0.001 1],'Type','real'); % [0.001 1] --> 0.5536
step_psi = optimizableVariable('step_psi',[1 20],'Type','integer'); % [0.1 10], real --> 9.9994
fun = @(x)runBehaviorDLDSSimStepPsiLambdaBehavior(dFF,behaviorData,...
    latentStatesX,groundTruthStates,simulatedPsi,x.lambda_behavior,x.step_psi);
results = bayesopt(fun,[lambda_behavior,step_psi],'Verbose',0,...
    'AcquisitionFunctionName','expected-improvement-plus');

% disp('Rewrite optimization based on mean of max varExplLatentsOnly for coeffs per ground truth coeff'...
%     'to answer the question of whether each gt coeff gets at least one good inferred coeff ')
% 
% nF = optimizableVariable('nF',[1 30],'Type','integer');
% lambda_behavior = optimizableVariable('lambda_behavior',[0.001 1],'Type','real');
% fun = @(x)runBehaviorDLDSSim_nFlambdabehavior(dFF,behaviorData,latentStatesX,x.nF,x.lambda_behavior);
% results = bayesopt(fun,[nF,lambda_behavior],'Verbose',0,...
%     'AcquisitionFunctionName','expected-improvement-plus');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Behavior and Psi - see after dLDS

%% Can coefficients be learned from behavior?

whichSample = 1;

for eachCoeff = 1:size(groundTruthStates{whichSample},1)
    predictorB = behaviorData{whichSample}.';
    outputOneC = groundTruthStates{whichSample}(eachCoeff,:).';
    [weights,fitinfo] = lasso(predictorB,outputOneC);
    [~,idxMSE] = min(fitinfo.MSE);
    figure()
    hold on
    plot(outputOneC)
    plot(weights(:,idxMSE).'*predictorB.')
    legend('actual c','predicted c')
end

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
inf_opts.lambda_historyb = 0.7;%0.45;           % ASC added 7/25 %0.7
inf_opts.tol             = 1e-8;           % 1e-3
inf_opts.max_iter2       = 20;             % 500 %20
inf_opts.max_iters       = input('max iters:');
inf_opts.F_update        = true; % default = true;
inf_opts.D_update        = false; % NoObs case - default = true;
inf_opts.N_ex            = 50;   
inf_opts.T_s             = 1000; 
inf_opts.step_d          = 1;  % ASC added 7/25
inf_opts.step_f          = 10;  % ASC added 7/25 % 30
inf_opts.step_decay      = 0.998;  
inf_opts.plot_option     = 10;  % ASC added 7/25
inf_opts.lambda_f_decay  = 0.996;
inf_opts.lambda_f        = 0.20;
inf_opts.solver_type     = ''; % fista
inf_opts.special         = 'noobs';
inf_opts.deltaDynamics   = false; %default: false %use x_t-x_{t-1}
inf_opts.AcrossIndividuals     = 0; %input('Across individuals? 1 for yes, 0 for no:');

inf_opts.behaviordLDS    = input('behavior_dLDS?:');
inf_opts.step_psi        = 1; % 30 for opt for psi
inf_opts.lambda_behavior = 0.5521; % 0.5536 for opt for psi


% inf_opts.solver_type     = 'tfocs'; %default: ''; (fista)
% inf_opts.CVX_Precision   = 'default';   %default = 'default'
% inf_opts.special         = ''; % regular bilinear inference
% inf_opts.debias          = true; % default = true; not used for NoObs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% [Phi, F] = bpdndf_dynamics_learning(dFF, [], simulatedDmatrix, inf_opts);              % Run the dictionary learning algorithm
% % Phi = simulatedDmatrix; F = Fgt2;
% [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
%                                        @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%%
if inf_opts.behaviordLDS
    [Phi, F, Psi] = bpdndf_dynamics_learning_behavior(dFF, [], [], [], behaviorData, inf_opts);              % Run the dictionary learning algorithm
else
    [Phi, F] = bpdndf_dynamics_learning(dFF, [], [], inf_opts);              % Run the dictionary learning algorithm
end

if inf_opts.behaviordLDS
    if inf_opts.AcrossIndividuals
        for ii = 1:size(Phi,1)
            [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference(dFF(ii), Phi{ii}, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
        end
    else
        [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
    
    end
else
    if inf_opts.AcrossIndividuals
        for ii = 1:size(Phi,1)
            [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference(dFF(ii), Phi{ii}, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
        end
    else
        [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
    
    end
end
%%

%%

if inf_opts.behaviordLDS
    plotMultiNetworkOutput(5, dFF, A_cell, B_cell, groundTruthStates, F, Psi, simulatedPsi)
else
    plotMultiNetworkOutput(5, dFF, A_cell, B_cell, groundTruthStates, F)
end

%% Can coefficients be learned from behavior?

whichSample = 1;

for eachCoeff = 1:size(B_cell{whichSample},1)
    predictorB = behaviorData{whichSample}.';
    outputOneC = B_cell{whichSample}(eachCoeff,:).';
    [weights,fitinfo] = lasso(predictorB,outputOneC);
    [~,idxMSE] = min(fitinfo.MSE);
    figure()
    hold on
    plot(outputOneC)
    plot(weights(:,idxMSE).'*predictorB.')
    legend('actual c','predicted c')
end

%%
varExplX = varExplLatentsOnly(A_cell,latentStatesX);
varExplX_summary = mean(varExplX);
%%
varExplC = varExplLatentsOnly(B_cell,groundTruthStates);

% Since there are different numbers of ground truth and inferred coeffs,
% each combination gets evaluated and returned for each sample. Then, we
% evaluate the maximum variance explained for each ground truth coeff among
% the inferred coeffs. Then we take the mean of these across all of the
% states and samples (in a cell array).

varExplC_summary = mean(cellfun(@(x) mean(max(x,[],1),"all","omitmissing"),varExplC), "all"); % should be close to 1


%% compare Psi to ground truth Psi

varExplPsi = varExplLatentsOnly({simulatedPsi.'},{Psi.'});
varExplPsi_summary = mean(cellfun(@(x) mean(max(x,[],1),"all","omitmissing"),varExplPsi), "all"); % should be close to 1
%%
figure()
set(gcf,"Color","white")
imagesc(simulatedPsi)
colorbar
title('Ground-truth Psi')

%%
figure()
set(gcf,"Color","white")
imagesc(Psi)
colorbar
title('Learned Psi')

%%
figure()
set(gcf,"Color","white")
subplot(2,4,1:4)
if exist('Psi','var')
    cmapMaxDistinct = distinguishable_colors(size(Psi,1));
else
    cmapMaxDistinct = distinguishable_colors(size(simulatedPsi,1));
end
% cmapMaxDistinct = distinguishable_colors(size(Psi,1));
hold on
tVals = 1:3000;
for ii=1:size(behaviorData{1},1)
    Y = behaviorData{1};
    plot(tVals,Y(ii,:),'color',cmapMaxDistinct(ii,:),'DisplayName',num2str(ii));
    hold on;
end
legend
ylabel('True behavior')
hold off

subplot(2,4,5:8)
hold on
for ii=1:size(Psi,1)
    Y = (Psi*B_cell{1});
    plot(tVals,Y(ii,:),'color',cmapMaxDistinct(ii,:),'DisplayName',num2str(ii));
    hold on;
end
legend
ylabel('Reconstructed behavior')
xlabel('Time')
hold off

%%
figure();
varExplPerBehaviorTrace = zeros(size(behaviorData,1),size(behaviorData{1},1));
for i = 1:size(behaviorData,1) % samples
    for j = 1:size(behaviorData{1},1) % behavior dimensions
        behaviorIJ = behaviorData{i}(j,:); % one behavior trace by time, from one sample
        reconstrBehaviorI = Psi*B_cell{i}; % all behavior traces by time, from one sample
        reconstrBehaviorIJ = reconstrBehaviorI(j,:);
        % % find spots where 1 frame increases, next frame decreases, or vice versa (unnatural jump)
        % idxIncrease = diff(reconstrBehaviorIJ)>1*std(reconstrBehaviorIJ);
        % idxDecrease = diff(reconstrBehaviorIJ)<-1*std(reconstrBehaviorIJ);
        % 
        % behaviorIJ(find(idxIncrease)) = behaviorIJ(find(idxIncrease)+1);
        % reconstrBehaviorIJ(find(idxIncrease)) = reconstrBehaviorIJ(find(idxIncrease)+1);
        % behaviorIJ(find(idxDecrease)) = behaviorIJ(find(idxDecrease)+1);
        % reconstrBehaviorIJ(find(idxDecrease)) = reconstrBehaviorIJ(find(idxDecrease)+1);
        % 
        % hold on
        % plot(behaviorIJ.')
        % plot(reconstrBehaviorIJ.')
        % legend('True behavior','Reconstructed behavior')
        % hold off
        % pause(1)
        % cla
        r = corrcoef(behaviorIJ,reconstrBehaviorIJ);
        varExplPerBehaviorTrace(i,j) = (r(1,2))^2;
    end
end

varExplBhv = mean(varExplPerBehaviorTrace,'all');
%%
figure()
set(gcf,"Color","white")
histogram(varExplPerBehaviorTrace,50)
title('Var Expl per behavior trace, all samples')

%%
r = corrcoef(behaviorData{1}.',(Psi*B_cell{1}).');
varExplBhv = (r(1,2))^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%