%% SSIM calculations 

addpath(genpath('/home/adam/Desktop/Versioned/CodeBaseRepo/trunk/CW_SSIM/'))
addpath(genpath('/home/adam/Desktop/Versioned/CodeBaseRepo/trunk/matlabPyrTools/'))

N_height = 4;
N_orient = 6;
guardb = 0;

for kk = 1:T_s
    %[vid_recon_cs(:), vid_recon_rwcs(:), vid_recon_dcs(:), vid_recon_drw(:), vid_recon_amp2(:), vid_recon_mdrw(:), vid_recon_modcs(:), vid_recon_rwl2b(:)];
    vid_ssim_cs(kk) = cssim_index_multi(255*(vid_recon_cs(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
    vid_ssim_rwcs(kk) = cssim_index_multi(255*(vid_recon_rwcs(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
    vid_ssim_dcs(kk) = cssim_index_multi(255*(vid_recon_dcs(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
    vid_ssim_drw(kk) = cssim_index_multi(255*(vid_recon_drw(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
    vid_ssim_amp2(kk) = cssim_index_multi(255*(vid_recon_amp2(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
    vid_ssim_mdrw(kk) = cssim_index_multi(255*(vid_recon_mdrw(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
    vid_ssim_modcs(kk) = cssim_index_multi(255*(vid_recon_modcs(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
    vid_ssim_mdrw(kk) = cssim_index_multi(255*(vid_recon_mdrw(:, :, kk) + 1)/2, ...
        255*(vid2(:, :, kk) + 1)/2, N_height, N_orient, guardb);
end

% for kk = 1:T_s
%     %[vid_recon_cs(:), vid_recon_rwcs(:), vid_recon_dcs(:), vid_recon_drw(:), vid_recon_amp2(:), vid_recon_mdrw(:), vid_recon_modcs(:), vid_recon_rwl2b(:)];
%     vid_ssim_cs(kk) = ssim_index(vid_recon_cs(:, :, kk), vid2(:, :, kk));
%     vid_ssim_rwcs(kk) = ssim_index(vid_recon_rwcs(:, :, kk), vid2(:, :, kk));
%     vid_ssim_dcs(kk) = ssim_index(vid_recon_dcs(:, :, kk), vid2(:, :, kk));
%     vid_ssim_drw(kk) = ssim_index(vid_recon_drw(:, :, kk), vid2(:, :, kk));
%     vid_ssim_amp2(kk) = ssim_index(vid_recon_amp2(:, :, kk), vid2(:, :, kk));
%     vid_ssim_mdrw(kk) = ssim_index(vid_recon_mdrw(:, :, kk), vid2(:, :, kk));
%     vid_ssim_modcs(kk) = ssim_index(vid_recon_modcs(:, :, kk), vid2(:, :, kk));
%     vid_ssim_rwl2b(kk) = ssim_index(vid_recon_rwl2b(:, :, kk), vid2(:, :, kk));
% end
%%
SSIM_ALL = [vid_ssim_cs(:), vid_ssim_rwcs(:), vid_ssim_dcs(:), vid_ssim_drw(:), vid_ssim_amp2(:), vid_ssim_mdrw(:), vid_ssim_modcs(:), vid_ssim_mdrw(:)];

% Time plot
line_style_list = {':b', '-g', '-.r', '--c', ':k', '--m', '-y', '--m'};
line_style_list2 = {'.b', '.g', '.r', '.c', '.k', '^m', 'oy', '^m'};
color_style_list = {'b', 'g', 'r', 'c', 'k', 'm', 'y', 'm'};

figure, hold on;
for kk = 1:8
    plot(mean(SSIM_ALL(:, kk), 3), line_style_list{kk}, 'LineWidth', 3)
%     plot(rMSE_ALL(:, kk, 2), line_style_list{kk}, 'LineWidth', 3)
end
box on
legend('Independent BPDN', 'Independent rw-BPDN', 'Dynamic BPDN', 'Dynamic rw-BPDN', 'DCS-AMP', 'RWL1-SS', 'modCS', 'WL1P')
set(gca, 'FontSize', 18, 'Xlim', [1,size(rMSE_ALL, 1)], 'Ylim', [0, 1])
xlabel('Frame Number', 'FontSize', 22)
ylabel('rMSE', 'FontSize', 22)
