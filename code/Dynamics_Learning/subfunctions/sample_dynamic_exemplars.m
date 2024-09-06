function x_im = sample_dynamic_exemplars(data_obj, opts)

% function x_im = sample_exemplars(data_obj, opts)
%
% Function to sample exemplars. Takes in a data object and a struct of
% algorithmic parameters and returns either a uniform sampling of the data
% or a sampling based on finding a sample set with low correlations between
% the samples. 
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

Ts    = opts.T_s;                                                          % Number of times
x_im  = cell(1,opts.N_ex);                                                 % Initialize the exemplar sample array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Training Data. 

if isnumeric(data_obj)&&(numel(size(data_obj))==2)

    data_use_ind = ceil((size(data_obj, 2)-Ts+1)*rand(1, opts.N_ex));      % Select a sampling of the data blindly
    
    for kk = 1:opts.N_ex
        x_im{kk} = data_obj(:,data_use_ind(kk):(data_use_ind(kk)+Ts-1));     % Extract a portion of the sequence
        if (opts.ssim_flag)&&(sqrt(sum(x_im{kk}(:).^2))>0)            % If desired, normalize the sequence
                x_im{kk} = x_im{kk}./sqrt(sum(x_im{kk}(:).^2));            % Normalize the chosen sample, as long as the norm is not zero
        end
    end
elseif isnumeric(data_obj)&&(numel(size(data_obj))==3)                     % Can't slice the main image cube, so extract data serially:   
    height_im = size(data_obj, 1) - opts.bl_size(1) + 1;                   % Get heights and widths for allowable starting indices for subimages
    width_im  = size(data_obj, 2) - opts.bl_size(2) + 1;                   %  ---------------
    num_im    = size(data_obj, 3) - Ts + 1;                                %  ---------------
    start_ind = ceil(rand(opts.N_ex, 3).*[height_im, width_im, num_im]);   % Pick random starting points
    end_ind   = bsxfun(@plus, start_ind, ...
                               [opts.bl_size(1), opts.bl_size(1), Ts]) - 1;% Calculate ending points
    for kk = 1:opts.N_ex
        x_im{kk} = reshape(data_obj(start_ind(kk,1):end_ind(kk,1),...
            start_ind(kk,2):end_ind(kk,2), start_ind(kk,3):end_ind(kk,3)), [], 1); % Pick out subimage and reshape it to a vector
        if (opts.ssim_flag)&&(sqrt(sum(x_im{kk}(:).^2))>0)            % If desired, normalize the sequence
                x_im{kk} = x_im{kk}./sqrt(sum(x_im{kk}(:).^2));            % Normalize the chosen sample, as long as the norm is not zero
        end
    end
elseif iscell(data_obj)&&(numel(size(data_obj{1}))==2)
    num_im    = length(data_obj);                                          % Get sizes of the data blocks (color images or video):
    samplesCounter = 0;
    timePointsSum = 0;
    for jj = 1:num_im
        timePointsSum = timePointsSum + size(data_obj{jj},2);
    end
    for jj = 1:num_im
        if opts.sampleProportionally
            trialProportion = size(data_obj{jj},2)/timePointsSum;
        else
            trialProportion = 1/num_im;
        end
        if jj == num_im
            nSamplesThisTrial = opts.N_ex - samplesCounter;
        else
            nSamplesThisTrial = ceil(trialProportion*opts.N_ex);
        end
        
        height_im = size(data_obj{jj}, 2) - Ts + 1;                             % Get heights and widths for allowable starting indices for subimages
        start_ind = ceil(bsxfun(@times,rand(nSamplesThisTrial,1),height_im));  % Pick random starting points
        for kk =  1:nSamplesThisTrial %opts.N_ex
            x_im{samplesCounter +kk} = data_obj{jj}(:,...
                                start_ind(kk):start_ind(kk)+Ts-1); % Pick out subimage and reshape it to a vector
            if (opts.ssim_flag)&&(sqrt(sum(x_im{samplesCounter +kk}(:).^2))>0)            % If desired, normalize the sequence
                    x_im{samplesCounter +kk} = x_im{samplesCounter +kk}./sqrt(sum(x_im{samplesCounter +kk}(:).^2));            % Normalize the chosen sample, as long as the norm is not zero
            end
        end
        samplesCounter = samplesCounter + nSamplesThisTrial;
    end
elseif iscell(data_obj)&&(numel(size(data_obj{1}))==3)
    disp('Note: sample_dynamic_examplars still uses dimensions of the first object in cell array for 3d images - must fix if you use this functionality')
    num_im = length(data_obj);                                             % Get sizes of the data blocks (color images or video):
    height_im = size(data_obj{1}, 1) - opts.bl_size(1) + 1;                % Get heights and widths for allowable starting indices for subimages
    width_im = size(data_obj{1}, 2) - opts.bl_size(2) + 1;
    depth_im = size(data_obj{1}, 3) - Ts + 1;
    
    start_ind = ceil(bsxfun(@times,rand(opts.N_ex, 4), ...
                                [height_im, width_im, depth_im, num_im])); % Pick random starting points
    end_ind   = bsxfun(@plus, start_ind(:,1:3), ...
                              [opts.bl_size(1), opts.bl_size(1), Ts]) - 1; % Calculate ending points
    for kk = 1:opts.N_ex
        x_im = reshape(data_obj{start_ind(kk,4)}(...
                                start_ind(kk,1):end_ind(kk,1),...
                                start_ind(kk,2):end_ind(kk,2), ...
                                start_ind(kk,3):end_ind(kk,3)), [], 1);    % Pick out subimage and reshape it to a vector
    end
else
    error('Unknown Data Type!! Choose ''vector'', ''square'' or ''cube''...')
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
