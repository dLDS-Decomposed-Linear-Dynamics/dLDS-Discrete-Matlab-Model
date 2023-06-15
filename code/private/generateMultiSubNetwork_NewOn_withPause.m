%% multi-subnetwork test case
% output: matrix of activity (time rows by channels
% cols)
% (Do we need ground truth?)
% What should signals look like? Should they just be sine waves of a
% certain frequency with noise? 
% How should they stack on top of each other? What do we mean by
% overlapping subnetworks?
% By a subnetwork, do we mean a set of channels that are active together?

function [xMatrix, varargout] = generateMultiSubNetwork_NewOn_withPause(nTp, nG1, nG2, decayOn, varargin)

nTimepoints = nTp;

nNeuronsGroup1 = nG1;
nNeuronsGroup2 = nG2;
nNeurons = nNeuronsGroup1 + nNeuronsGroup2;

xMatrix = ones(nNeurons,nTimepoints);% initialize matrix nTimepoints by nNeurons
% prevent the values from exploding
% find the eigenvalues of the random matrix
% standardize the eigenvalues and reconstruct U * lambda * Uinv

if ~isempty(varargin)
    f1Stdzd = varargin{1};
    f2Stdzd = varargin{2};
    f3Stdzd = varargin{3};
else
    f1 = randn(nNeuronsGroup1, nNeuronsGroup1);
    [U1,L1] = eig(f1);
    f1Stdzd = real(U1 * (L1.*diag(1./max(abs(L1))))/U1); % pay attention: inverse, not transpose
    f2 = randn(nNeuronsGroup2, nNeuronsGroup2);
    [U2,L2] = eig(f2);
    f2Stdzd = real(U2 * (L2.*diag(1./max(abs(L2))))/U2);
    f3 = randn(nNeuronsGroup2, nNeuronsGroup2);
    [U3,L3] = eig(f3);
    f3Stdzd = real(U3 * (L3.*diag(1./max(abs(L3))))/U3);
end

if decayOn
    t = (1:1:nTimepoints)';     % seconds
    primeVec = primes(100);
    indicesPrimes = randsample(1:size(primeVec,2),3);
    freq2 = 1/primeVec(indicesPrimes(2));
    freq3 = 1/primeVec(indicesPrimes(3));
    f2coeffs = 0.1*sin(2*pi*freq2*t)+1;
    f3coeffs = 0.1*sin(2*pi*freq3*t)+1;
    %concatenate
    coeffsTwoOnDecay = [ones(nTimepoints,nNeuronsGroup1) repmat(f2coeffs,[1,nNeuronsGroup2])];
    coeffsNewOnDecay = [ones(nTimepoints,nNeuronsGroup1) repmat(f3coeffs,[1,nNeuronsGroup2])];
end


topRows = cat(2,f1Stdzd, zeros(nNeuronsGroup1,nNeuronsGroup2));
bottomRows1 = cat(2,zeros(nNeuronsGroup2, nNeuronsGroup1), f2Stdzd);
bottomRows2 = cat(2,zeros(nNeuronsGroup2, nNeuronsGroup1), f3Stdzd);
bottomRowsDefault = zeros(nNeuronsGroup2,nNeurons);


fMatrixOneOn = [topRows; bottomRowsDefault];
fMatrixTwoOn = [topRows; bottomRows1];
fMatrixNewOn = [topRows; bottomRows2];

for i = 1:round(nTimepoints/5)-1
    xMatrix(:,i+1) = fMatrixOneOn * xMatrix(:,i);
end
for i = round(nTimepoints/5):round(2*nTimepoints/5)-1
    if i == round(nTimepoints/5)
        fakeLastStep = [xMatrix(1:nNeuronsGroup1,i); ones(nNeuronsGroup2,1)];
        xMatrix(:,i+1) = fMatrixTwoOn * fakeLastStep;
    else
        xMatrix(:,i+1) = fMatrixTwoOn * xMatrix(:,i);
    end
end
for i = round(2*nTimepoints/5):round(3*nTimepoints/5)-1
    if decayOn
        xMatrix(:,i+1) = diag(coeffsTwoOnDecay(i,:))*fMatrixTwoOn * xMatrix(:,i);
    else
        xMatrix(:,i+1) = fMatrixOneOn * xMatrix(:,i);
    end
end
for i = round(3*nTimepoints/5):round(4*nTimepoints/5)-1
    if i == round(3*nTimepoints/5)
        fakeLastStep = [xMatrix(1:nNeuronsGroup1,i); ones(nNeuronsGroup2,1)];
        xMatrix(:,i+1) = fMatrixNewOn * fakeLastStep;
    else
        xMatrix(:,i+1) = fMatrixNewOn * xMatrix(:,i);
    end
end
for i = round(4*nTimepoints/5):nTimepoints
    if decayOn
        xMatrix(:,i+1) = diag(coeffsTwoOnDecay(i,:))*fMatrixNewOn * xMatrix(:,i);
    else
        xMatrix(:,i+1) = fMatrixOneOn * xMatrix(:,i);
    end
end

if nargout > 1
    varargout{1} = {f1Stdzd, f2Stdzd, f3Stdzd};
    if decayOn
        varargout{2} = {f2coeffs, f3coeffs};
    end
end
end