%% multi-subnetwork test case
% output: matrix of activity (time rows by channels
% cols)
% (Do we need ground truth?)
% What should signals look like? Should they just be sine waves of a
% certain frequency with noise? 
% How should they stack on top of each other? What do we mean by
% overlapping subnetworks?
% By a subnetwork, do we mean a set of channels that are active together?

function [xMatrix, varargout] = generateMultiSubNetwork(nTp, nG1, nG2)

nTimepoints    = nTp;
nNeuronsGroup1 = nG1;
nNeuronsGroup2 = nG2;
nNeurons       = nNeuronsGroup1 + nNeuronsGroup2;

xMatrix = ones(nNeurons,nTimepoints);% initialize matrix nTimepoints by nNeurons
% prevent the values from exploding
% find the eigenvalues of the random matrix
% standardize the eigenvalues and reconstruct U * lambda * Uinv
f1 = randn(nNeuronsGroup1, nNeuronsGroup1);
[U1,L1] = eig(f1);
f1Stdzd = real(U1 * (L1.*diag(1./max(abs(L1))))/U1); % pay attention: inverse, not transpose
f2 = randn(nNeuronsGroup2, nNeuronsGroup2);
[U2,L2] = eig(f2);
f2Stdzd = real(U2 * (L2.*diag(1./max(abs(L2))))/U2);
topRows = cat(2,f1Stdzd, zeros(nNeuronsGroup1,nNeuronsGroup2));
bottomRows = cat(2,zeros(nNeuronsGroup2, nNeuronsGroup1), f2Stdzd);
fMatrixAllOn = [topRows; bottomRows];
fMatrixOneGroupOff = [topRows; zeros(nNeuronsGroup2, nNeurons)];

for i = 1:round(nTimepoints/2)-1
    xMatrix(:,i+1) = fMatrixAllOn * xMatrix(:,i);
end

for i = round(nTimepoints/2):nTimepoints-1
    xMatrix(:,i+1) = fMatrixOneGroupOff * xMatrix(:,i);
end

if nargout > 1
    varargout{1} = {f1Stdzd, f2Stdzd};
end
end
