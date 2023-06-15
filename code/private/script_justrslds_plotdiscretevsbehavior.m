%%
addpath(genpath('../../..'))

%% load data
% load('saveCElegansWorkspace_v4_afterrestoreddynamicslearning.mat');
dat = load('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegansData/WT_Stim.mat');

rsldsStimWorm1_10x4z = load('Stim_Worm1_discretestates_10x4z.mat');
rsldsStimWorm2_10x4z = load('Stim_Worm2_discretestates_10x4z.mat');
rsldsStimWorm3_10x4z = load('Stim_Worm3_discretestates_10x4z.mat');
rsldsStimWorm4_10x4z = load('Stim_Worm4_discretestates_10x4z.mat');
rsldsStimWorm5_10x4z = load('Stim_Worm5_discretestates_10x4z.mat');
rsldsStimWorm6_10x4z = load('Stim_Worm6_discretestates_10x4z.mat');
rsldsStimWorm7_10x4z = load('Stim_Worm7_discretestates_10x4z.mat');

rsldsStimWorms_10x4z(1,:) = {rsldsStimWorm1_10x4z};
rsldsStimWorms_10x4z(2,:) = {rsldsStimWorm2_10x4z};
rsldsStimWorms_10x4z(3,:) = {rsldsStimWorm3_10x4z};
rsldsStimWorms_10x4z(4,:) = {rsldsStimWorm4_10x4z};
rsldsStimWorms_10x4z(5,:) = {rsldsStimWorm5_10x4z};
rsldsStimWorms_10x4z(6,:) = {rsldsStimWorm6_10x4z};
rsldsStimWorms_10x4z(7,:) = {rsldsStimWorm7_10x4z};

rsldsStimWorm1_10x10z = load('Stim_Worm1_discretestates_10x10z.mat');
rsldsStimWorm2_10x10z = load('Stim_Worm2_discretestates_10x10z.mat');
rsldsStimWorm3_10x10z = load('Stim_Worm3_discretestates_10x10z.mat');
rsldsStimWorm4_10x10z = load('Stim_Worm4_discretestates_10x10z.mat');
rsldsStimWorm5_10x10z = load('Stim_Worm5_discretestates_10x10z.mat');
rsldsStimWorm6_10x10z = load('Stim_Worm6_discretestates_10x10z.mat');
rsldsStimWorm7_10x10z = load('Stim_Worm7_discretestates_10x10z.mat');

rsldsStimWorms_10x10z(1,:) = {rsldsStimWorm1_10x10z};
rsldsStimWorms_10x10z(2,:) = {rsldsStimWorm2_10x10z};
rsldsStimWorms_10x10z(3,:) = {rsldsStimWorm3_10x10z};
rsldsStimWorms_10x10z(4,:) = {rsldsStimWorm4_10x10z};
rsldsStimWorms_10x10z(5,:) = {rsldsStimWorm5_10x10z};
rsldsStimWorms_10x10z(6,:) = {rsldsStimWorm6_10x10z};
rsldsStimWorms_10x10z(7,:) = {rsldsStimWorm7_10x10z};


%% reorder rslds classifications to keep neighbor states close
% wormdiscrete = rsldsStimWorms_10x4z(1);
% wormdiscrete = wormdiscrete{1,1}.discretestates;
% neighborMat = areTheyNeighbors(wormdiscrete);
    
%% plot rslds and behavior (7 stim worms)

totalBhvStates = 4;

figure();

for i = 1:7
    dFF = dat.WT_Stim(i).traces;
    bhv = dat.WT_Stim(i).States;
    tVals    = dat.WT_Stim(i).timeVectorSeconds;
    frames = length(tVals);
    wormdiscrete = rsldsStimWorms_10x4z(i);
    wormdiscrete = wormdiscrete{1,1}.discretestates;
%     wormdiscretesparse = rsldsStimWorms_10x4z_sparse(i,:,:);
    wormdiscrete = reorderNeighbors(wormdiscrete);
    subplot(7,2,(i-1)*2+1)
    hold on

    bottom = 0;
    top = totalBhvStates+1;
    sz = 8 .* ones(frames,1);
    yyaxis left
%     scatter(tVals, bhv, sz, 'MarkerFaceAlpha',0.5);
    scatter(tVals, bhv, sz, '.');
    ylabel('Behavioral state');
    ylim('padded')
    
    yyaxis right
%     scatter(tVals, wormdiscrete, sz, 'x');
    plot(tVals, wormdiscrete);
    ylabel('rSLDS label');
    legend('Behavior','rSLDS');
    xlabel('Time (s)');
    ylim('padded')    
    set(gca, 'TickDir', 'out');
    box off;
    hold off;

end

for i = 1:7
    dFF = dat.WT_Stim(i).traces;
    bhv = dat.WT_Stim(i).States;
    tVals    = dat.WT_Stim(i).timeVectorSeconds;
    frames = length(tVals);
    wormdiscrete = rsldsStimWorms_10x10z(i);
    wormdiscrete = wormdiscrete{1,1}.discretestates;
%     wormdiscretesparse = rsldsStimWorms_10x10z_sparse(i,:,:);
    wormdiscrete = reorderNeighbors(wormdiscrete);
    subplot(7,2,i*2)
    hold on

    bottom = 0;
    top = totalBhvStates+1;
    sz = 8 .* ones(frames,1);
    yyaxis left
%     scatter(tVals, bhv, sz, 'MarkerFaceAlpha',0.5);
    scatter(tVals, bhv, sz, '.');
    ylabel('Behavioral state');
    ylim('padded')

    
    yyaxis right
%     scatter(tVals, wormdiscrete, sz, 'x');
    plot(tVals, wormdiscrete);
    ylabel('rSLDS label');
    legend('Behavior','rSLDS');
    xlabel('Time (s)');
    ylim('padded')
    set(gca, 'TickDir', 'out');
    box off;
    hold off;

end

%% 
function reorderedData = reorderNeighbors(labelsToReorder)
% This function returns the rslds discrete states vector with the labels
% renumbered so that each state that usually follows another is ordered
% accordingly.

doNotReorder = false;

if doNotReorder
    reorderedData = labelsToReorder;
else
    nOriginalStates = max(labelsToReorder)+1;
%     neighborlinessMatrix = zeros(nOriginalStates,nOriginalStates);
    reorderedData = zeros(size(labelsToReorder));
    deltaT = 100; % Set time window size for looking ahead to which states follow.
    neighborlinessMatrix = createNeighborlinesMatrix(labelsToReorder,nOriginalStates,deltaT);

    
    % Reweight columns (second values in pairs) by base rates.
    for j = 1:nOriginalStates
        baseRateOfLabelJ = sum(labelsToReorder(:)==(j-1))/size(labelsToReorder,2); % because rslds labels start from 0
        neighborlinessMatrix(:,j) = round((1/baseRateOfLabelJ) * neighborlinessMatrix(:,j));
    end

    
    [bestRows,bestCols] = findIndsBestPairs(neighborlinessMatrix,nOriginalStates);
    reorderedStateLabels = reorderTheStateLabels(nOriginalStates,bestRows,bestCols);
   
    % Return the sorted data.
    % Remember to add 1 to the rslds label value because rslds labels start
    % from 0 (python). Use the original rslds label (+1) as the index to 
    % pick the new state label from the reordered vector.
    for i = 1:size(labelsToReorder,2)
        reorderedData(i) = reorderedStateLabels(labelsToReorder(i)+1);
    end
end
end

function neighborsMat = createNeighborlinesMatrix(originalLabels, nStates, timeWindowSize)
neighborsMat = zeros(nStates,nStates);
for i = 1:nStates
    for k = 1:nStates
        for t = 1:size(originalLabels,2)-timeWindowSize
            if originalLabels(t) == i-1 % because rslds labels start from 0
                window = t:t+timeWindowSize; % asymmetric % Set time window size for looking ahead to which states follow.
                neighborsMat(i,k) = neighborsMat(i,k)...
                    + sum(originalLabels(1,window)==(k-1)); % because rslds labels start from 0
            end
        end
    end
end
end

function [indBestPairsRows, indBestPairsCols] = findIndsBestPairs(neighborsMat,nStates)
% Find the row and column indices of the n-1 best pairs in the matrix
% (without redundancy).
% Set diagonal of neighborsMat to 0 because we don't need
% self-pairing scores.
neighborsMat = neighborsMat - diag(diag(neighborsMat));
% Initialize
indBestPairsRows = zeros(1,nStates-1);
indBestPairsCols = zeros(1,nStates-1);
for i = 1:nStates-1
    [~,indMaxVal] = max(neighborsMat(:));  %index of top value in whole matrix (down cols, then over rows)
    % Turn index back into row,col format.
    modulusIndStates = mod(indMaxVal,nStates);
    if modulusIndStates == 0
        indRow = nStates;
    else
        indRow = modulusIndStates;
    end
    indMaxValDouble = double(indMaxVal);
    nOriginalStatesDouble = double(nStates);
    indCol = ceil(double(indMaxValDouble)/double(nOriginalStatesDouble));
    
    % Store the indices of this best value in row, col order
    indBestPairsRows(i) = indRow;
    indBestPairsCols(i) = indCol;

    % If we set 1 (indRow) followed by 2 (indCol) with the highest
    % neighborliness, we won't want to later set 1 followed by another 
    % value. Neither do we want to set another value followed by 2.
    % Therefore, we remove that row and column from the matrix.
    neighborsMat(indRow,:) = zeros(1,nStates);
    neighborsMat(indCol,:) = zeros(nStates,1);
    
end
end

function newStateLabels = reorderTheStateLabels(nStates,indBestPairsRows,indBestPairsCols)
    % Permute order of columns in asymmetric neighborlinessMatrix.
    % original order [1,2,3,4]
    newStateLabels = (1:nStates);
    for i = nStates-1:-1:1
        % Set pairs so indices next to each other, and set the next column's
        % index to where the paired column was.
        % Sort by worst of the top pairings first, 
        % to prioritize the best pairings (overriding the weaker ones).
        firstLabelInPair = indBestPairsRows(i);
        secondLabelInPair = indBestPairsCols(i);
        currentIndexOfFirstLabelInPair ...
            = find(newStateLabels == firstLabelInPair);
        currentIndexOfSecondLabelInPair ...
            = find(newStateLabels == secondLabelInPair);
        % Accounting for the size of the label vector:
        % if you want to move a label to the end, shift back all the
        % values to fill the gap in the middle first.
        if currentIndexOfFirstLabelInPair ~= nStates
            % Copy over the old label that followed your first label in the
            % pair of interest to the spot where the second label in the 
            % pair is; then put the second label in the pair right after 
            % the first label in the pair.
            newStateLabels(currentIndexOfSecondLabelInPair) ...
                = newStateLabels(currentIndexOfFirstLabelInPair+1);
            newStateLabels(currentIndexOfFirstLabelInPair+1) ...
                = secondLabelInPair;
        elseif currentIndexOfFirstLabelInPair == nStates
            % Shift all the values after the old spot of the
            % second label in the pair back by the difference between the 
            % last index in the vector and the old index of the second  
            % label in the pair.
            shiftBackBy = currentIndexOfFirstLabelInPair-currentIndexOfSecondLabelInPair; 
            newStateLabels(currentIndexOfSecondLabelInPair:currentIndexOfSecondLabelInPair+shiftBackBy-1) ...
                = newStateLabels(currentIndexOfSecondLabelInPair+1:currentIndexOfSecondLabelInPair+shiftBackBy); 
            newStateLabels(nStates) = secondLabelInPair;
        else
            disp('How did you get here?')
        end
        disp(newStateLabels)
    end
end

% function bestPermutation = findBestPermutation(nStates,neighborsMat)
% %idea: instead of sorting, just calculate total neighborliness score for
% %each permutation and pick highest
% permsStates = perms(1:nStates);
% totalNeighborlinessOfEachPerm = zeros(size(permsStates,1),1);
% for i = 1:size(permsStates,1)
%     for j = 1:nStates-1
%         pairStatesScore = neighborsMat(permsStates(i,j),permsStates(i,j+1)); %find neighborliness value from the matrix for the pair presented in the permutation
%         totalNeighborlinessOfEachPerm(i) = totalNeighborlinessOfEachPerm(i) + pairStatesScore; %add that pair's score to the total for that permutation
%     end
% end
% 
% [~,indexMaxNeighborliness] = max(totalNeighborlinessOfEachPerm);
% 
% bestPermutation = permsStates(indexMaxNeighborliness(1),:);
% 
% end