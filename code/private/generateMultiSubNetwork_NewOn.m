function [xMatrix, varargout] = generateMultiSubNetwork_NewOn(nTp, nG1, nG2, nTr, varargin)

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
nTrials        = nTr;

nNeuronsGroup1 = nG1;
nNeuronsGroup2 = nG2;
nNeurons       = nNeuronsGroup1 + nNeuronsGroup2;

if nargin > 4;   OneFOnWholeTrial = varargin{1};
else;            OneFOnWholeTrial = false;
end

if nargin > 5;   TwoFOnWholeTrial = varargin{2};
else             TwoFOnWholeTrial = false;    
end

if nargin > 6;   Multiswitch = varargin{3};
else             Multiswitch = false;    
end

% prevent the values from exploding
% find the eigenvalues of the random matrix
% standardize the eigenvalues and reconstruct U * lambda * Uinv
f1Stdzd = sampRandDyn(nNeuronsGroup1);
f2Stdzd = sampRandDyn(nNeuronsGroup2);
f3Stdzd = sampRandDyn(nNeuronsGroup2);
f4Stdzd = sampRandDyn(nNeuronsGroup1);

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
    xMatrix = runTwoSystemSimulation(nNeurons,nTimepoints,nTr,topRows,bottomRows1,bottomRows2);
    if nargout > 1
        varargout{1} = {f1Stdzd, f2Stdzd, f3Stdzd};
    end
elseif Multiswitch
    xMatrix = cell(nTr,1);
    for ll = 1:nTr
        xMatrix{ll} = ones(nNeurons,nTimepoints);
        xMatrix{ll}(:,1) = randn(nNeurons,1);
        tNow = 0;
        while tNow < nTimepoints-1
            tJump = round(rand(1)*300)+100;
            tJump = tJump - max(tNow + tJump-nTimepoints+1,0);
            fMatrix = mkFmat(topRows, topRows2, bottomRows1, bottomRows2);
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
    if nargout > 1
        varargout{1} = {f1Stdzd, f2Stdzd, f3Stdzd, f4Stdzd};
    end
end


end


function f1Stdzd = sampRandDyn(nNeuronsGroup)

f = randn(nNeuronsGroup, nNeuronsGroup);
[U,L] = eig(f);
f1Stdzd = real(U*(L.*diag(1./max(abs(L))))/U); % pay attention: inverse, not transpose

end

function fMatrix = mkFmat(topRows, topRows2, bottomRows1, bottomRows2)

topOpt = round(2*rand(1));
botOpt = round(2*rand(1));

fMatrix = [(topOpt==1)*topRows     + (topOpt==0)*topRows2; 
           (botOpt==1)*bottomRows1 + (botOpt==0)*bottomRows2];

end

function xMatrix = runTwoSystemSimulation(nNeurons,nTimepoints,nTr,topRows,bottomRows1,bottomRows2)

    fMatrixAllOn = [topRows; bottomRows1];
    fMatrixNewOn = [topRows; bottomRows2];

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