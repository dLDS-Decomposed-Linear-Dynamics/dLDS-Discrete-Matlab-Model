%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
% lambdaCombos = [[1,1];[1,2];[2,1];[2,2];[2,4];[4,2];[4,4];[1,8];[8,1];[8,8]];
% lambdaCombos = [[1,1];[1,2];[2,2];[2,4];[4,4];[1,8];[8,8]];
lambdaCombos = [1,1];
nCombos = size(lambdaCombos,1);
% fstruct = dir('*_lambdaval_*_lambdab_*.mat');

workspacesLambdaCombos = cell(1, nCombos);

for k = 1:nCombos
  myfilename = sprintf('saveMSNWorkspace_3000i_10i2_Finitsvddecay_Kp5_Lval_p%d_Lb_p%d_tol_1em3_D10.mat', lambdaCombos(k,1), lambdaCombos(k,2));
  workspacesLambdaCombos{k} = importdata(myfilename);
end

rsldsWorm1S = load('./Worm1_WT_Stim_Zhat_discretestates.mat');

%% Precomputation of useful things


% tVals    = workspacesLambdaCombos{1}.dat.WT_Stim(1).timeVectorSeconds;
% tVals    = workspacesLambdaCombos{1}.dat2.WT_NoStim(5).timVectorSeconds;

% c_values = B_cell{1}.';
% x_values = A_cell{1};
% x_values_timebyx  = x_values.';
% dataReconstructed = (Phi * x_values).';
% 
% frames = size(x_values,2);

ifHighlight = false;

if ifHighlight

    wormdiscrete = rsldsWorm1S.discretestatesWorm1S;
    wormdiscretesparse = zeros(size(wormdiscrete,1),1);
    for rowI = 1:size(wormdiscrete,1)
        wormdiscretesparse(rowI) = find(wormdiscrete(rowI,:)==max(wormdiscrete(rowI,:)),1);
    end

    % flag points (indices) of oscillation/uncertainty in rslds
    flaggedSuspectStateClassification = zeros(size(wormdiscretesparse,1),1);
    for i = 2:size(wormdiscretesparse,1)-1
        if (wormdiscretesparse(i) ~= wormdiscretesparse(i-1)) && (wormdiscretesparse(i-1) == wormdiscretesparse(i+1))
            flaggedSuspectStateClassification(i) = 1;
        end
    end

    verticesFlaggedOn = zeros(size(wormdiscretesparse,1),1);
    verticesFlaggedOff = zeros(size(wormdiscretesparse,1),1);

    % define patches of oscillation in terms of their vertices (start and end indices)

    for j = 11:size(flaggedSuspectStateClassification,1)-10
        if (flaggedSuspectStateClassification(j) == 1) && (sum(flaggedSuspectStateClassification(j+1:j+10))>2) && (sum(flaggedSuspectStateClassification(j-10:j-1))==0)
            verticesFlaggedOn(j) = 1;
        elseif (flaggedSuspectStateClassification(j) == 1) && (sum(flaggedSuspectStateClassification(j-10:j-1))>2) && (sum(flaggedSuspectStateClassification(j+1:j+10))==0)
            verticesFlaggedOff(j) = 1;
        end
    end

    % correction - close off the last patch
    if size(find(verticesFlaggedOn==1)) > size(find(verticesFlaggedOff==1))
        verticesFlaggedOff(end) = 1;
    end

    idxVerticesOn  = find(verticesFlaggedOn==1);
    idxVerticesOff = find(verticesFlaggedOff==1);
    idxVerticesOn  = idxVerticesOn*max(tVals)./numel(tVals);
    idxVerticesOff = idxVerticesOff*max(tVals)./numel(tVals);


    % common patch params
    grey  = [127 127 127]./255;
    colorF = grey; % face color
    alphaF = 0.1; % face opacity

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot dynamics and states

figure();
for k = 1:nCombos
    subplot(nCombos,1,k);

    c_values = workspacesLambdaCombos{k}.B_cell{1}.';
    x_values = workspacesLambdaCombos{k}.A_cell{1};
    x_values_timebyx  = x_values.';
    thisDFF = workspacesLambdaCombos{k}.dFF;
    dataReconstructed = (workspacesLambdaCombos{k}.Phi * x_values).'; % note: not rescaled
    
    frames = size(x_values,2);
    
    disp('tVals = frames - for simulated');
    tVals = 1:frames; 

    hold on;
    bottom = min(c_values,[],'all');
    top = max(c_values,[],'all');
    if ifHighlight
        for i = 1:size(idxVerticesOn)
            vertices = [idxVerticesOn(i) bottom; idxVerticesOff(i) bottom; idxVerticesOff(i) top; idxVerticesOn(i) top];
            faces = [1     2     3     4];
            patch('Faces',faces,'Vertices',vertices,'FaceColor',colorF, 'FaceAlpha', alphaF,'EdgeColor','none');
        end
    end
    plot(tVals, c_values); 
    xlabel('Time (s)');
    axis([0 max(tVals) bottom top]);
    title(sprintf('Estimated coefficient c values, lambda val 0.%d lambda b 0.%d', lambdaCombos(k,1), lambdaCombos(k,2)));
    set(gca, 'TickDir', 'out');
    box off;
    hold off;

%     % note: not rescaled
%     flatDFF = reshape(thisDFF, numel(thisDFF),[]);
%     flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
%     rval = corrcoef(flatDFF, flatReconstructed);
%     varExpl = (rval(1,2))^2;
%     disp(varExpl);
end

%% 
figure();
for k = 1:nCombos
    subplot(nCombos,1,k);

    c_values = workspacesLambdaCombos{k}.B_cell{1}.';
    x_values = workspacesLambdaCombos{k}.A_cell{1};
    x_values_timebyx  = x_values.';
    thisDFF = workspacesLambdaCombos{k}.dFF;
    dataReconstructed = (workspacesLambdaCombos{k}.Phi * x_values).'; % note: not rescaled
    
    frames = size(x_values,2);
    
    disp('tVals = frames - for simulated');
    tVals = 1:frames; 

    hold on;
    bottom = min(x_values,[],'all');
    top = max(x_values,[],'all');
    if ifHighlight
        for i = 1:size(idxVerticesOn)
            vertices = [idxVerticesOn(i) bottom; idxVerticesOff(i) bottom; idxVerticesOff(i) top; idxVerticesOn(i) top];
            faces = [1     2     3     4];
            patch('Faces',faces,'Vertices',vertices,'FaceColor',colorF, 'FaceAlpha', alphaF,'EdgeColor','none');
        end
    end
    plot(tVals, x_values); 
    xlabel('Time (s)');
    axis([0 max(tVals) bottom top]);
    title(sprintf('Estimated latent states x, lambda val 0.%d lambda b 0.%d', lambdaCombos(k,1), lambdaCombos(k,2)));
    set(gca, 'TickDir', 'out');
    box off;
    hold off;

    % note: not rescaled
    flatDFF = reshape(thisDFF, numel(thisDFF),[]);
    flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
    rval = corrcoef(flatDFF, flatReconstructed);
    varExpl = (rval(1,2))^2;
    disp(varExpl);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




