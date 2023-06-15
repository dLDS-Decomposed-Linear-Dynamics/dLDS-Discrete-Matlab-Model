function [xMatrix, varargout] = generateContMultiSubNetwork(nTp, nG, nSC, nTr, varargin)

%% multi-subnetwork test case
% output: matrix of activity (time rows by channels
% cols)
% (Do we need ground truth?)
% What should signals look like? Should they just be sine waves of a
% certain frequency with noise? 
% How should they stack on top of each other? What do we mean by
% overlapping subnetworks?
% By a subnetwork, do we mean a set of channels that are active together?

nTimepoints    = nTp;
nNeuronsGroup1 = nG(1);
nNeuronsGroup2 = nG(2);
nSubCirc1      = nSC(1);
nSubCirc2      = nSC(2);
nNeurons       = nNeuronsGroup1 + nNeuronsGroup2;

if nargin > 4;   contOpt = varargin{1};
else;            contOpt = false;    
end

if nargin > 5;   OneFOnWholeTrial = varargin{2};
else;            OneFOnWholeTrial = false;
end

if nargin > 6;   TwoFOnWholeTrial = varargin{3};
else;            TwoFOnWholeTrial = false;    
end

if nargin > 7;   Multiswitch = varargin{4};
else;            Multiswitch = false;    
end

% prevent the values from exploding
% find the eigenvalues of the random matrix
% standardize the eigenvalues and reconstruct U * lambda * Uinv
fSubCirc1 = cell(nSubCirc1,1);
fSubCirc2 = cell(nSubCirc2,1);

for ll = 1:nSubCirc1
    fSubCirc1{ll} = sampRandDyn(nNeuronsGroup1, contOpt);
end

for ll = 1:nSubCirc2
    fSubCirc2{ll} = sampRandDyn(nNeuronsGroup2, contOpt);
    if OneFOnWholeTrial
        fSubCirc2{ll} = 0*fSubCirc2{ll}; 
    elseif TwoFOnWholeTrial
        if ll > 1; fSubCirc2{ll} = 0*fSubCirc2{ll}; end
    end
end

topRows    = cell(nSubCirc1,1);
bottomRows = cell(nSubCirc2,1);
for ll = 1:nSubCirc1
    topRows{ll}    = cat(2,fSubCirc1{ll}, zeros(nNeuronsGroup1,nNeuronsGroup2));
end
for ll = 1:nSubCirc2
    bottomRows{ll} = cat(2,zeros(nNeuronsGroup2, nNeuronsGroup1), fSubCirc2{ll});
end


if ~Multiswitch
    xMatrix = runTwoSystemSimulation(nNeurons,nTimepoints,nTr,topRows,bottomRows, contOpt);
    if nargout > 1
        varargout{1} = {fSubCirc1, fSubCirc2};
    end
elseif Multiswitch
    if contOpt
        [xMatrix, groudTruthC] = generateRandSwitchSysCont(nTr,nNeurons,...
                               nTimepoints, topRows, bottomRows, contOpt);
    else
        [xMatrix, groudTruthC] = generateRandSwitchSysDisc(nTr, ...
                         nNeurons, nNeuronsGroup1, nNeuronsGroup2, ...
                               nTimepoints, topRows, bottomRows, contOpt);
    end
    if nargout > 1;  varargout{1} = {fSubCirc1, fSubCirc2};     end
    if nargout > 2;  varargout{2} = groudTruthC;     end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f1Stdzd = sampRandDyn(nNeuronsGroup, varargin)

if nargin > 1; contOpt = varargin{1};
else;          contOpt = true;            end

f1Stdzd = RandOrthMat(nNeuronsGroup, 1e-5);

if contOpt
    [U,L]   = eig(f1Stdzd);
    f1Stdzd = real(U*diag(log(diag(L)))/U);                           % pay attention: inverse, not transpose
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fMatrix, varargout] = mkFmat(topRows, bottomRows, contOpt)

nTop = numel(topRows);
nBot = numel(bottomRows);

topOpt = floor((1+nTop)*rand(1));
botOpt = floor((1+nBot)*rand(1));

if contOpt; dirVals = 2*(rand(1,2)>0.5)-1;
else;       dirVals = [1,1];
end

if topOpt > 0;   topSel = dirVals(1)*topRows{topOpt};
else;            topSel = 0*topRows{1};
end
if botOpt > 0;   botSel = dirVals(2)*bottomRows{botOpt};
else;            botSel = 0*bottomRows{1};
end

fMatrix = [topSel; botSel];
if nargout > 1; varargout{1} = [topOpt*dirVals(1), botOpt*dirVals(2)]; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xMatrix = runTwoSystemSimulation(nNeurons,nTimepoints,nTr,topRows,bottomRows, contOpt)

    if contOpt
        fMatrixAllOn = expm([topRows{1}; bottomRows{1}]);
        fMatrixNewOn = expm([topRows{1}; bottomRows{2}]);
    end
    
    if nTr == 1
        xMatrix = ones(nNeurons,nTimepoints);
        for i = 1:round(nTimepoints/2)-1
            xMatrix(:,i+1) = fMatrixAllOn * xMatrix(:,i);
        end
        for i = round(nTimepoints/2):nTimepoints-1
            xMatrix(:,i+1) = fMatrixNewOn * xMatrix(:,i);
        end
    else
        xMatrix = cell(nTr,1);
        for ll = 1:nTr
            xMatrix{ll}      = ones(nNeurons,nTimepoints);
            xMatrix{ll}(:,1) = randn(nNeurons,1);
            for i = 1:round(nTimepoints/2)-1
                xMatrix{ll}(:,i+1) = fMatrixAllOn*xMatrix{ll}(:,i);
            end
            for i = round(nTimepoints/2):nTimepoints-1
                xMatrix{ll}(:,i+1) = fMatrixNewOn*xMatrix{ll}(:,i);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xMatrix, groudTruthC] = generateRandSwitchSysDisc(nTr, nNeurons, nNeuronsGroup1, nNeuronsGroup2, nTimepoints, topRows, bottomRows, contOpt)

groudTruthC = cell(nTr,1);
xMatrix = cell(nTr,1);
for ll = 1:nTr
    xMatrix{ll}      = ones(nNeurons,nTimepoints);
    xMatrix{ll}(:,1) = randn(nNeurons,1);
    groudTruthC{ll}  = zeros(numel(topRows)+numel(bottomRows),nTimepoints);
    tNow             = 0;
    while tNow < nTimepoints-1
        tJump = round(rand(1)*300)+100;
        tJump = tJump - max(tNow + tJump-nTimepoints+1,0);
        [fMatrix, selF] = mkFmat(topRows, bottomRows, contOpt);
        if all(xMatrix{ll}(1:nNeuronsGroup1,tNow+1)==0)
            xMatrix{ll}(1:nNeuronsGroup1,tNow+1)=randn(nNeuronsGroup1,1);
        end
        if all(xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)==0)
            xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)=randn(nNeuronsGroup2,1);
        end
        
        xMatrix{ll}(1:nNeuronsGroup1,tNow+1)=xMatrix{ll}(1:nNeuronsGroup1,tNow+1)./norm(xMatrix{ll}(1:nNeuronsGroup1,tNow+1));
        xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)=xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)./norm(xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1));
        
        for i = (tNow+1):(tNow+tJump)
            xMatrix{ll}(:,i+1) = fMatrix*xMatrix{ll}(:,i);
        end
        
        if selF(1) > 0
            groudTruthC{ll}(selF(1), (tNow+1):(tNow+tJump+1)) = 1;
        end
        if selF(2) > 0
            groudTruthC{ll}(numel(topRows)+selF(2), (tNow+1):(tNow+tJump+1)) = 1;
        end
        
        tNow  = tNow + tJump;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xMatrix, groudTruthC] = generateRandSwitchSysCont(nTr, nNeurons, nTimepoints, topRows, bottomRows, contOpt)

groudTruthC = cell(nTr,1);
xMatrix     = cell(nTr,1);

for ll = 1:nTr
    xMatrix{ll}      = ones(nNeurons,nTimepoints);
    xMatrix{ll}(:,1) = randn(nNeurons,1);
    xMatrix{ll}(:,1) = xMatrix{ll}(:,1)./norm(xMatrix{ll}(:,1));
    groudTruthC{ll}  = zeros(numel(topRows)+numel(bottomRows),nTimepoints);
    
    tNow             = 0;
    while tNow < nTimepoints-1
        tJump   = round(rand(1)*300)+100;
        tJump   = tJump - max(tNow + tJump-nTimepoints+1,0);
        [fMatrix, selF] = mkFmat(topRows, bottomRows, contOpt);
        fMatrix = expm(fMatrix);
        for i = (tNow+1):(tNow+tJump+1)
            xMatrix{ll}(:,i+1) = fMatrix*xMatrix{ll}(:,i);
        end
        if selF(1) > 0
            groudTruthC{ll}(selF(1), (tNow+1):(tNow+tJump+1)) = 1;
        end
        if selF(1) > 0
            groudTruthC{ll}(selF(1), (tNow+1):(tNow+tJump+1)) = 1;
        end
        tNow  = tNow + tJump;
        
        if selF(1) > 0
            groudTruthC{ll}(abs(selF(1)), (tNow+1):(tNow+tJump+1)) = sign(selF(1));
        end
        if selF(2) > 0
            groudTruthC{ll}(numel(topRows)+abs(selF(2)), (tNow+1):(tNow+tJump+1)) = sign(selF(2));
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%