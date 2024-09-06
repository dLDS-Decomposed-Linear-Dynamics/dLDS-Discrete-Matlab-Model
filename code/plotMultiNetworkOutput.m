function plotMultiNetworkOutput(tSel, dFF, A_cell, B_cell, groundTruthStates, F,varargin)
if nargin > 6
    nLines = 7;
    Psi = varargin{1};
    simulatedPsi = varargin{2};
else
    nLines = 5;
end


% tSel = 5; % selected trial
figure(1), cla; 
subplot(nLines,4,[1,4]), plot(dFF{tSel}.')       % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('Signal (GT)')
set(gca,'XLim',[1,size(dFF{tSel},2)]);
subplot(nLines,4,[5,8]), plot(A_cell{tSel}.')    % ASC added 7/25
set(gca,'XLim',[1,size(dFF{tSel},2)]);
box off; xlabel('Time (frames)'); ylabel('Reconstruction')

subplot(nLines,4,[9,12]), imagesc(abs(B_cell{tSel}))   % ASC added 7/25
set(gca,'XLim',[1,size(dFF{tSel},2)]);
box off; xlabel('Time (frames)'); ylabel('Coeff Index')

% subplot(5,4,[9,12]), plot(B_cell{tSel}.')   % ASC added 7/25
% set(gca,'XLim',[1,size(dFF{tSel},2)]);
% box off; xlabel('Time (frames)'); ylabel('Dynamics coefficient amplitude')
% legend('C1','C2','C3','C4')
subplot(nLines,4,[13,16]), imagesc(groundTruthStates{tSel})   % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('True state'); ylim('padded')
%set(gca,'XLim',[1,size(dFF{tSel},2)],'YLim',[-0.1,1.1]);
% legend('C1','C2','C3','C4')

allFs = F{1};
nF1   = size(F{1},1);
for ll = 2:numel(F)
    allFs = cat(2,allFs,[zeros(nF1,5),F{ll}]);
end
subplot(nLines,4,[17,20]), imagesc(allFs); 
axis image; axis off; colormap gray

if nargin > 6
    subplot(nLines,4,[21,24]), imagesc(Psi); 
    subplot(nLines,4,[25,28]), imagesc(simulatedPsi); 
end

% subplot(5,4,17), imagesc(F{1}); axis image; axis off; colormap gray
% subplot(5,4,18), imagesc(F{2}); axis image; axis off; colormap gray
% subplot(5,4,19), imagesc(F{3}); axis image; axis off; colormap gray
% subplot(5,4,20), imagesc(F{4}); axis image; axis off; colormap gray

end