%% multi-subnetwork test case 2
% output: matrix of activity (time rows by channels cols)
% with constant overlapping activation between 2 groups of neurons

function [xMatrix, varargout] = generateMultiSubNetwork_Overlap(nTp, nG1, nG2)

nTimepoints = nTp;

nNeuronsGroup1 = nG1;
nNeuronsGroup2 = nG2;
nNeurons = nNeuronsGroup1 + nNeuronsGroup2;

xMatrix = ones(nNeurons,nTimepoints);% initialize matrix nTimepoints by nNeurons
% prevent the values from exploding
% find the eigenvalues of the random matrix
% standardize the eigenvalues and reconstruct U * lambda * Uinv
f1 = randn(nNeuronsGroup1, nNeuronsGroup1);
% [U1,L1] = eig(f1);
% f1Stdzd = real(U1 * (L1.*diag(1./max(abs(L1))))/U1); % pay attention: inverse, not transpose
f2 = randn(nNeuronsGroup2, nNeuronsGroup2);
% [U2,L2] = eig(f2);
% f2Stdzd = real(U2 * (L2.*diag(1./max(abs(L2))))/U2);
topRows = cat(2,f1, zeros(nNeuronsGroup1,nNeuronsGroup2));
bottomRows = cat(2,zeros(nNeuronsGroup2, nNeuronsGroup1), f2);
fMatrixAllOn = [topRows; bottomRows];
% fMatrixOneGroupOff = [topRows; zeros(nNeuronsGroup2, nNeurons)];

f3 = randn(5,5);

fMatrixAllOnOverlap = fMatrixAllOn;
fMatrixAllOnOverlap(nNeuronsGroup1-2:nNeuronsGroup1+2, nNeuronsGroup1-2:nNeuronsGroup1+2) = f3;

[U3,L3] = eig(fMatrixAllOnOverlap);
f3Stdzd = real(U3 * (L3.*diag(1./max(abs(L3))))/U3);

fMatrixAllOnOverlap = f3Stdzd;
fMatrixAllOnOverlap(nNeuronsGroup1+1:nNeuronsGroup1+2,1:nNeuronsGroup1-3) = 0;
fMatrixAllOnOverlap(1:nNeuronsGroup1-3,nNeuronsGroup1+1:nNeuronsGroup1+2) = 0;
fMatrixAllOnOverlap(nNeuronsGroup1+3:end, 1:nNeuronsGroup1) = 0;
fMatrixAllOnOverlap(1:nNeuronsGroup1, nNeuronsGroup1+3:end) = 0;

fMatrixOneGroupOffOverlap = fMatrixAllOnOverlap;
fMatrixOneGroupOffOverlap(nNeuronsGroup1+3:end, :) = 0;
fMatrixOneGroupOffOverlap(:,nNeuronsGroup1+3:end) = 0;

    
for i = 1:round(nTimepoints/2)-1
    xMatrix(:,i+1) = fMatrixAllOnOverlap * xMatrix(:,i);
end

for i = round(nTimepoints/2):nTimepoints-1
    xMatrix(:,i+1) = fMatrixOneGroupOffOverlap * xMatrix(:,i);
end

figure();
imagesc(xMatrix);
colorbar;

1+1;


% if nargout > 1
%     varargout{1} = {f1Stdzd, f2Stdzd};
% end
end