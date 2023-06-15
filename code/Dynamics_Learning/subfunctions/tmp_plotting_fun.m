function tmp_plotting_fun(data_type,D,F,inf_opts,n_iters,dataShape,varargin)

% tmp_plotting_fun(data_type,D,F,inf_opts,varargin)
%
% Function to plot the learned representational and dynamics dictionaries.
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 6
    D_true = varargin{1};
else
    D_true = zeros(size(D));
end

if nargin > 7
    F_true = varargin{2};
else
    F_true = cell(size(F));
    for kk = 1:numel(F)
        F_true{kk} = zeros(size(F{kk}));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform plotting

if strcmp(data_type, 'synth')
    if mod(n_iters, inf_opts.plot_option) == 0
        figure(1)
        subplot(2,2,1), imagesc(basis2img2(D, sqrt(sig_opts.N)*[1,1], sqrt(sig_opts.N)*[1,1]))
        title('Learned Dictionary', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        subplot(2,2,2), imagesc((D)*F{1}*(D'))
        title('Learned Dynamics', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        figure(1)
        subplot(2,2,3), imagesc(basis2img2(D_true, sqrt(sig_opts.N)*[1,1], sqrt(sig_opts.N)*[1,1]))
        title('True Dictionary', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        subplot(2,2,4), imagesc(F_true{1})
        title('True Dynamics', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        drawnow
    end
elseif strcmp(data_type, 'BBC')||strcmp(data_type,'datamatrix')
    Fmat = zeros(numel(F{1}), numel(F));
    if mod(n_iters, inf_opts.plot_option) == 0
        %%
        figure(1)
        if strcmp(dataShape, 'image')
            subplot(1,2,1), imagesc(basis2img2(D, sqrt(inf_opts.M)*[1,1], sqrt(inf_opts.N)*[1,1]))
            title('Learned Dictionary', 'FontSize', 20)
            colormap gray
            axis image
            axis off
        else
        end
        for mm = 1:numel(F)
                Fmat(:, mm) = reshape(F{mm}, [], 1);
        end
        subplot(1,2,2), imagesc(basis2img2(Fmat, size(F{1}), [5,5]))
        title('Learned Dynamics', 'FontSize', 20)
        colormap gray
        axis image
        axis off
        drawnow
    end
else
end


end