vid_blur = zeros(size(vid2).*[0.5, 0.5, 1]);

for kk = 1:size(vid_blur, 3)
    vid_blur(:, :, kk) = blur_operator(vid2(:, :, kk)) + sqrt(noise_var)*randn(sqrt(M), sqrt(M));
end
%%

OUT_VID = VideoWriter('Video_Test.avi');
OUT_VID.FrameRate = 20; %,'Uncompressed AVI');
for nn = 1:size(vid2, 3)
    figure(1)
    subplot(2, 4, 1), imagesc(vid_blur(:, :, nn)) %vid2 %imagesc(reshape(x_store(:, nn), 24, 24))
    title('True Signal', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    subplot(2, 4, 2), imagesc(vid_recon_cs(:, :, nn))
    title('BPDN', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    subplot(2, 4, 3), imagesc(vid_recon_dcs(:, :, nn))
    title('BPDN-DF', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    subplot(2, 4, 4), imagesc(vid_recon_drw(:, :, nn))
    title('RWL1-DF', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    drawnow
    subplot(2, 4, 5), imagesc(vid_recon_rwcs(:, :, nn))
    title('RWL1', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    subplot(2, 4, 6), imagesc(vid_recon_mdrw(:, :, nn))
    title('modRWL1', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    subplot(2, 4, 7), imagesc(vid_recon_modcs(:, :, nn))
    title('modCS', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    subplot(2, 4, 8), imagesc(vid_recon_amp2(:, :, nn))
    title('DCS-AMP', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    drawnow
    set(gcf, 'color', [1,1,1])
    OUT_MOV(nn) = getframe(gcf);
%         pause(0.5)
end
open(OUT_VID);
writeVideo(OUT_VID, OUT_MOV);
close(OUT_VID);
