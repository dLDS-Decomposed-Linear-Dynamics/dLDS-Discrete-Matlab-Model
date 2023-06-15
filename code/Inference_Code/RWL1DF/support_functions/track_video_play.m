function track_video_play(x, varargin)

if nargin > 1
    x_size = varargin{1};
else
    x_size = sqrt(size(x, 1))*ones(1, 2);
end

if nargin > 2
    fig_no = varargin{2};
else
    fig_no = [];
end

figure(fig_no)
for kk = 1:size(x, 2)
    imagesc(reshape(x(:, kk), x_size))
    colormap gray
    axis off
    axis image
    pause(0.1)
end

end