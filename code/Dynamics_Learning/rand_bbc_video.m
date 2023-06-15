function X_ex = rand_bbc_video(sel_dims, varargin)

% X_ex = rand_bbc_video(sel_dims, varargin)
% 
% Sample BBC videos randomly
% 
% 
% 2015 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing and Error checking

if nargin == 1
    downsamp_factor = 1;
    chunk_path = '/home/adam/GITrepos/dynamics_learning/data/cadieu_bbc_video/BBCmotion/orig/chunk';
else
    downsamp_factor = varargin{1};
    chunk_path = varargin{2};
end

if numel(sel_dims) == 1
    sel_dims = sel_dims*[1,1,1];
elseif numel(sel_dims) == 2
    sel_dims = [sel_dims(1),sel_dims(1),sel_dims(2)];
else                                                                       % Nothing to do
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for down-sampling 

if downsamp_factor == 1                                                    % no downsampling
else
%     sel_dims = [downsamp_factor*sel_dims(1),downsamp_factor*sel_dims(2),sel_dims(3)];
    sel_dims = [sel_dims(1),sel_dims(2),sel_dims(3)];                      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Randomly select and load video
vid_select = ceil(80*rand(1));                                             % Select a random chuck to load
X = [];
load(sprintf('%s%d', chunk_path, vid_select),'X')

dim_start(3) = ceil(rand(1)*(size(X,3) - sel_dims(3)+1));

X = double(X(:,:,dim_start(3):(dim_start(3)+sel_dims(3)-1)));
Xd = X(1:downsamp_factor:end, 1:downsamp_factor:end, 1);

X_tmp2 = zeros(size(Xd,1),size(Xd,2), size(X,3));                          % Initialize new array that is the right size
if downsamp_factor == 1                                                    % no downsampling
else
    for kk = 1:size(X, 3)
        X_tmp2(:,:,kk) = imresize(X(:,:,kk),1/downsamp_factor,'lanczos3'); % Perform downsampling
    end
end

vid_size = size(X_tmp2);                                                   % Make sure the right video size is used


% Extract random portion of the video
dim_start(1) = ceil(rand(1)*(vid_size(1) - sel_dims(1)+1));
dim_start(2) = ceil(rand(1)*(vid_size(2) - sel_dims(2)+1));


X_tmp = X(dim_start(1):(dim_start(1)+sel_dims(1)-1), dim_start(2):(dim_start(2)+sel_dims(2)-1) ,:);

X_ex = zeros(sel_dims(1)*sel_dims(2), sel_dims(3));
for kk = 1:size(X_tmp, 3)
    X_ex(:, kk) = reshape(X_tmp(:, :, kk) - mean(mean(double(X_tmp(:, :, kk)))), [], 1);
    if max(abs(X_ex(:, kk))) == 0
        fprintf('WARNING: all zero video found\n')
        X_ex(:, kk) = X_ex(:, kk);
    else
        X_ex(:, kk) = X_ex(:, kk)/max(abs(X_ex(:, kk)));
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%