%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Script to analyze BBC dynamics dictionary %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data

load ~/GITrepos/dynamics_learning/results/OLD/20150218_BBCvideo_12x12_4x_nF20.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize dictionary and dynamics

I = basis2img2(D, [12,12], 2*[12,12]);
figure(1);
imagesc(I);
colormap gray
axis image
axis off
set(gcf, 'color', [1,1,1])

figure(2);
subplot(1,2,1), imagesc(F{2});
colormap gray
axis image
axis off
subplot(1,2,2), imagesc(F{3});
colormap gray
axis image
axis off
set(gcf, 'color', [1,1,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

for kk = 1:numel(F)
    Fsums(:,kk) = sum(abs(F{kk}),1);
end

figure; imagesc(Fsums.')
% axis image
axis off
colormap gray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% C240, F1&4 OR F1&5 OR F1&6
% C3, F2&4 OR 20&4
% C495 F4&10 
% C396 F1&12
% C575 F1&23
% C574 F1&2

im = zeros(24);
im(end-4) = 1;
im = reshape(im, [],1);
Nmax = 10;
figure(3); set(gcf, 'color', [1,1,1]);

fsel = [8,5];
DT   = 1;

% subplot(1+Nmax/DT,5,1), imagesc(reshape(D*im, [12,12]));
% colormap gray
% axis image
% axis off
% title('Original Dictionary Element', 'FontSize', 18)    

for kk = 0:DT:Nmax
    for ll = 0:4
        a = ll/sqrt(ll^2 + (4-ll)^2);
        b = (4-ll)/sqrt(ll^2 + (4-ll)^2);
        im_now = reshape(D*((a*F{fsel(1)} + b*F{fsel(2)})^kk)*im, [12,12]);
        subplot(1+Nmax/DT,5,1+5*kk/DT + ll), imagesc(im_now);
        colormap gray
        axis image
        axis off
%         title(sprintf('Frame %d', kk), 'FontSize', 18)    
    end
    drawnow
    pause(0.1)
end



%%

% C572 F8&5
figure(4); set(gcf, 'color', [1,1,1]);
for ll = 0:19
    a = ll/sqrt(ll^2 + (19-ll)^2);
    b = (19-ll)/sqrt(ll^2 + (19-ll)^2);
    im_now = reshape(D*((a*F{fsel(1)} + b*F{fsel(2)})^1)*im, [12,12]);
    subplot(2,10,1 + ll), imagesc(im_now);
    colormap gray
    axis image
    axis off
%     title(sprintf('Frame %d', kk), 'FontSize', 18)    
end
drawnow
pause(0.1)


%%

figure(5); set(gcf, 'color', [1,1,1]);
Nmax = 20;
Ncol = 1+Nmax/DT;
for kk = 0:DT:Nmax
    for ll = 1:5
%         a = (1+ll)/sqrt((1+ll)^2 + (6-ll)^2);
%         b = (6-ll)/sqrt((1+ll)^2 + (6-ll)^2);
        a = (3 + ll)/12;
        b = (12 - 3 - ll)/12;
        im_now = reshape(D*((a*F{fsel(1)} + b*F{fsel(2)})^kk)*im, [12,12]);
        subplot(5,Ncol,1 + kk + (ll-1)*Ncol), cla;
        subplot(5,Ncol,1 + kk + (ll-1)*Ncol), imagesc(im_now/max(im_now(:)));
        colormap gray
        axis image
        axis off
%         title(sprintf('Frame %d', kk), 'FontSize', 18)    
    end
    drawnow
    pause(0.1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

U      = cell(size(F));
Lam    = cell(size(F));
f_rank = zeros(size(F));
for kk = 1:numel(F)
    [U{kk}, Lam{kk}] = eig(F{kk});
    Lam{kk}          = diag(Lam{kk});
    f_rank(kk)       = rank(F{kk});
end

%%

figure(6)
subplot(2,2,1), histogram2(real(reshape(cell2mat(Lam), [],1)), imag(reshape(cell2mat(Lam), [],1)))
set(gca, 'TickDir','out','FontSize',18,'XLim',[-1,1],'YLim',[-1,1])


U_tmp = cell2mat(U);
U_tmp = U_tmp(:,[1:576:end, 2:576:end]);
U_corr = abs(U_tmp'*U_tmp);
subplot(2,2,2), histogram(U_corr(triu(ones(size(U_corr)),1)==1))
box off
set(gca, 'TickDir','out','FontSize',18)
xlabel('Top two eigenvector correlations','FontSize', 18)
set(gcf, 'color', [1,1,1]);

LamSum = bsxfun(@times,cumsum(abs(cell2mat(Lam)).^2),1./sum(abs(cell2mat(Lam)).^2));
subplot(2,2,3), histogram(sum(LamSum(:,1:24)<0.99))
box off
% set(gca, 'Yscale', 'log', 'Xscale', 'log')
 
subplot(2,2,4), imagesc(F{end})
colormap gray
axis image
axis off

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
