function [xMatrix, varargout] = generateContMultiSubNetwork_storeStates(nTp, nG1, nG2, nTr, varargin)

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
nNeuronsGroup1 = nG1;
nNeuronsGroup2 = nG2;
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
f1Stdzd = sampRandDyn(nNeuronsGroup1, contOpt);
f2Stdzd = sampRandDyn(nNeuronsGroup2, contOpt);
f3Stdzd = sampRandDyn(nNeuronsGroup2, contOpt);
f4Stdzd = sampRandDyn(nNeuronsGroup1, contOpt);

if OneFOnWholeTrial
    f2Stdzd = 0*f2Stdzd;
    f3Stdzd = 0*f3Stdzd;
else
    if TwoFOnWholeTrial;   f3Stdzd = 0*f2Stdzd;    end
end

topRows     = cat(2,f1Stdzd, zeros(nNeuronsGroup1,nNeuronsGroup2));
topRows2    = cat(2,f4Stdzd, zeros(nNeuronsGroup1,nNeuronsGroup2));
bottomRows1 = cat(2,zeros(nNeuronsGroup2, nNeuronsGroup1), f2Stdzd);
bottomRows2 = cat(2,zeros(nNeuronsGroup2, nNeuronsGroup1), f3Stdzd);

if ~Multiswitch
    xMatrix = runTwoSystemSimulation(nNeurons,nTimepoints,nTr,topRows,bottomRows1,bottomRows2, contOpt);
    if nargout > 1
        varargout{1} = {f1Stdzd, f2Stdzd, f3Stdzd};
    end
elseif Multiswitch
    if contOpt
        [xMatrix, groundTruthStates] = generateRandSwitchSysCont(nTr, nNeurons, ...
                        nNeuronsGroup1, nNeuronsGroup2, nTimepoints, ...
                        topRows, topRows2, bottomRows1, bottomRows2);
    else
        [xMatrix, groundTruthStates] = generateRandSwitchSysDisc(nTr, nNeurons, ...
                      nNeuronsGroup1, nNeuronsGroup2, nTimepoints, ...
                           topRows, topRows2, bottomRows1, bottomRows2);
    end
    if nargout > 1
        varargout{1} = {f1Stdzd, f2Stdzd, f3Stdzd, f4Stdzd};
    end
    if nargout > 2
        varargout{2} = groundTruthStates;
    end
end


end


function f1Stdzd = sampRandDyn(nNeuronsGroup, varargin)

if nargin > 1; contOpt = varargin{1};
else;          contOpt = true;            end

f1Stdzd = RandOrthMat(nNeuronsGroup, 1e-5);

if contOpt
    [U,L]   = eig(f1Stdzd);
    f1Stdzd = real(U*diag(log(diag(L)))/U);                           % pay attention: inverse, not transpose
end

end

function [fMatrix, groundTruthStates] = mkFmat(topRows, topRows2, bottomRows1, bottomRows2, groundTruthStates, tNow, tJump, varargin)

if nargin > 7; contOpt = varargin{1};
else;          contOpt = false;            end

if contOpt
    topOpt = round(2*rand(1));
    botOpt = round(2*rand(1));

    fMatrix = [(topOpt==1)*topRows*(rand()<0.5)     + (topOpt==0)*topRows2*(rand()<0.5); 
               (botOpt==1)*bottomRows1*(rand()<0.5) + (botOpt==0)*bottomRows2*(rand()<0.5)];

    % FIXME: ground truth tracker
else
    topOpt = round(2*rand(1));
    botOpt = round(2*rand(1));

    fMatrix = [(topOpt==1)*topRows     + (topOpt==0)*topRows2; 
               (botOpt==1)*bottomRows1 + (botOpt==0)*bottomRows2];
    
    if topOpt == 1
        groundTruthStates(1,(tNow+1):(tNow+tJump)) = 1;
    elseif topOpt == 0
        groundTruthStates(4,(tNow+1):(tNow+tJump)) = 4;
    end %otherwise leave as 0

    if botOpt == 1
        groundTruthStates(2,(tNow+1):(tNow+tJump)) = 2;
    elseif botOpt == 0
        groundTruthStates(3,(tNow+1):(tNow+tJump)) = 3;
    end %otherwise leave as 0

end


end

function xMatrix = runTwoSystemSimulation(nNeurons,nTimepoints,nTr,topRows,bottomRows1,bottomRows2, contOpt)

    if contOpt
        fMatrixAllOn = expm([topRows; bottomRows1]);
        fMatrixNewOn = expm([topRows; bottomRows2]);
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

function [xMatrix,groundTruthStates] = generateRandSwitchSysDisc(nTr, nNeurons, nNeuronsGroup1, nNeuronsGroup2, nTimepoints, topRows, topRows2, bottomRows1, bottomRows2)

xMatrix = cell(nTr,1);
groundTruthStates = cell(nTr,1);
for ll = 1:nTr
    groundTruthStates{ll} = zeros(1,nTimepoints);

    xMatrix{ll}      = ones(nNeurons,nTimepoints);
    xMatrix{ll}(:,1) = randn(nNeurons,1);
    tNow             = 0;
    while tNow < nTimepoints-1
        tJump = round(rand(1)*300)+100;
        tJump = tJump - max(tNow + tJump-nTimepoints+1,0);
        [fMatrix,groundTruthStates{ll}] = mkFmat(topRows, topRows2, bottomRows1, bottomRows2, groundTruthStates{ll}, tNow,tJump);
        if all(xMatrix{ll}(1:nNeuronsGroup1,tNow+1)==0)
            xMatrix{ll}(1:nNeuronsGroup1,tNow+1)=randn(nNeuronsGroup1,1);
        end
        if all(xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)==0)
            xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)=randn(nNeuronsGroup2,1);
        end
        for i = (tNow+1):(tNow+tJump)
            xMatrix{ll}(:,i+1) = fMatrix*xMatrix{ll}(:,i);
        end
        tNow  = tNow + tJump;
    end
end
end


function [xMatrix,groundTruthStates] = generateRandSwitchSysCont(nTr, nNeurons, nNeuronsGroup1, nNeuronsGroup2, nTimepoints, topRows, topRows2, bottomRows1, bottomRows2)

xMatrix = cell(nTr,1);
for ll = 1:nTr
    groundTruthStates{ll} = zeros(1,nTimepoints);

    xMatrix{ll}      = ones(nNeurons,nTimepoints);
    xMatrix{ll}(:,1) = randn(nNeurons,1);
    xMatrix{ll}(:,1) = xMatrix{ll}(:,1)./norm(xMatrix{ll}(:,1));
    tNow             = 0;
    while tNow < nTimepoints-1
        tJump   = round(rand(1)*300)+100;
        tJump   = tJump - max(tNow + tJump-nTimepoints+1,0);
        [fMatrix,groundTruthStates{ll}] = mkFmat(topRows, topRows2, bottomRows1, bottomRows2, groundTruthStates{ll}, tNow,tJump);
        fMatrix = expm(fMatrix);
%         if all(xMatrix{ll}(1:nNeuronsGroup1,tNow+1)==0)
%             xMatrix{ll}(1:nNeuronsGroup1,tNow+1)=randn(nNeuronsGroup1,1);
%         end
%         if all(xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)==0)
%             xMatrix{ll}(nNeuronsGroup1+(1:nNeuronsGroup2),tNow+1)=randn(nNeuronsGroup2,1);
%         end
        for i = (tNow+1):(tNow+tJump+1)
            xMatrix{ll}(:,i+1) = fMatrix*xMatrix{ll}(:,i);
        end
        tNow  = tNow + tJump;
    end
end
end

