%% Track Code

% Calculate where targets are

thresh_val_cs = 0.24;
thresh_val_dcs = 0.139;
thresh_val_drw = 0.72;

x_loc = x_store ~= 0;
x_cs_loc = abs(x_cs_all) > thresh_val_cs;
x_out_loc = abs(x_out_all) > thresh_val_dcs;
x_lsm1_loc = abs(x_lsm1_all) > thresh_val_drw;


cs_image_total = zeros(24, 24, 3);
dcs_image_total = zeros(24, 24, 3);
drw_image_total = zeros(24, 24, 3);

OUT_VID = VideoWriter('Track_Test_Locate_single.avi');
OUT_VID.FrameRate = 5; %,'Uncompressed AVI');
for nn = 1:size(x_cs_all, 2)
    figure(1)
    subplot(1, 4, 1), imagesc(reshape(x_loc(:, nn), 24, 24))
    title('True Signal', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    cs_image = zeros(24, 24, 3);
    cs_image(:, :, 3) = reshape(x_cs_loc(:, nn)&x_loc(:, nn), 24, 24);
    cs_image(:, :, 2) = reshape((~x_cs_loc(:, nn))&x_loc(:, nn), 24, 24);
    cs_image(:, :, 1) = reshape(x_cs_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    cs_image_total = cs_image_total|cs_image;
    
    subplot(1, 4, 2), imagesc(cs_image_total) %reshape(x_cs_loc(:, nn), 24, 24))
    title('BPDN', 'FontSize', 25)
    xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_cs_loc(:, nn))), ...
        sum((~x_loc(:, nn))&x_cs_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    
    dcs_image = zeros(24, 24, 3);
    dcs_image(:, :, 3) = reshape(x_out_loc(:, nn)&x_loc(:, nn), 24, 24);
    dcs_image(:, :, 2) = reshape((~x_out_loc(:, nn))&x_loc(:, nn), 24, 24);
    dcs_image(:, :, 1) = reshape(x_out_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    dcs_image_total = dcs_image_total|dcs_image;
    
    subplot(1, 4, 3), imagesc(dcs_image_total) %reshape(x_out_loc(:, nn), 24, 24))
    title('BPDN-DF', 'FontSize', 25)
    xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_out_loc(:, nn))), ...
        sum((~x_loc(:, nn))&x_out_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    
    drw_image = zeros(24, 24, 3);
    drw_image(:, :, 3) = reshape(x_lsm1_loc(:, nn)&x_loc(:, nn), 24, 24);
    drw_image(:, :, 2) = reshape((~x_lsm1_loc(:, nn))&x_loc(:, nn), 24, 24);
    drw_image(:, :, 1) = reshape(x_lsm1_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    drw_image_total = drw_image_total|drw_image;
    
    subplot(1, 4, 4), imagesc(drw_image_total) %reshape(x_lsm1_loc(:, nn), 24, 24))
    title('RWL1-DF', 'FontSize', 25)
    xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_lsm1_loc(:, nn))), ...
        sum((~x_loc(:, nn))&x_lsm1_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    drawnow
    set(gcf, 'color', [1,1,1])
    OUT_MOV(nn) = getframe(gcf);
    % pause(0.5)
end
open(OUT_VID);
writeVideo(OUT_VID, OUT_MOV);
close(OUT_VID);

fprintf('Probability of Missed Detection: %f, %f, %f \n', ...
    sum(sum(x_loc&(~x_lsm1_loc), 1)), sum(sum(x_loc&(~x_cs_loc), 1)), ...
    sum(sum(x_loc&(~x_out_loc), 1)))
fprintf('Probability of False Alarm: %f, %f, %f \n', ...
    sum(sum((~x_loc)&x_lsm1_loc, 1)), sum(sum((~x_loc)&x_cs_loc, 1)), ...
    sum(sum((~x_loc)&x_out_loc, 1)))
    
%%



figure(2)
subplot(1, 2, 1), plot((sum(x_loc&(~x_lsm1_loc), 1)), '-b', 'LineWidth', 3)
subplot(1, 2, 1), hold on
plot((sum(x_loc&(~x_cs_loc), 1)), '--r', 'LineWidth', 3)
plot((sum(x_loc&(~x_out_loc), 1)), ':k', 'LineWidth', 3)
title('Missed Detections', 'FontSize', 25)
set(gca, 'FontSize', 20)
hold off

subplot(1, 2, 2), plot((sum((~x_loc)&x_lsm1_loc, 1)), '-b', 'LineWidth', 3)
subplot(1, 2, 2), hold on
plot((sum((~x_loc)&x_cs_loc, 1)), '--r', 'LineWidth', 3)
plot((sum((~x_loc)&x_out_loc, 1)), ':k', 'LineWidth', 3)
title('False Alarms', 'FontSize', 25)
set(gca, 'FontSize', 20)
hold off

%%

figure(2)
subplot(1, 2, 1), plot(cumsum(sum(x_loc&(~x_lsm1_loc), 1)), '-b', 'LineWidth', 3)
subplot(1, 2, 1), hold on
plot(cumsum(sum(x_loc&(~x_cs_loc), 1)), '--r', 'LineWidth', 3)
plot(cumsum(sum(x_loc&(~x_out_loc), 1)), '.-k', 'LineWidth', 3)
title('Missed Detections', 'FontSize', 25)
set(gca, 'FontSize', 20)
hold off

subplot(1, 2, 2), plot(cumsum(sum((~x_loc)&x_lsm1_loc, 1)), '-b', 'LineWidth', 3)
subplot(1, 2, 2), hold on
plot(cumsum(sum((~x_loc)&x_cs_loc, 1)), '--r', 'LineWidth', 3)
plot(cumsum(sum((~x_loc)&x_out_loc, 1)), '.-k', 'LineWidth', 3)
title('False Alarms', 'FontSize', 25)
set(gca, 'FontSize', 20)
hold off

%%

color_limits = [min([min(min(min(x_store))), min(min(min(x_store - x_cs_all))), ...
    min(min(min(x_store - x_out_all))), min(min(min(x_store - x_lsm1_all)))]), ...
    max([max(max(max(x_store))), max(max(max(x_store - x_cs_all))), ...
    max(max(max(x_store - x_out_all))), max(max(max(x_store - x_lsm1_all)))])];

OUT_VID = VideoWriter('Track_Test_Error_noRMSE.avi');
OUT_VID.FrameRate = 5; %,'Uncompressed AVI');
for nn = 1:size(x_cs_all, 2)
    figure(1)
    subplot(1, 4, 1), imagesc(reshape(x_store(:, nn), 24, 24))
    title('True Signal', 'FontSize', 25)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)
    
    subplot(1, 4, 2), imagesc(reshape((x_store(:, nn) - x_cs_all(:, nn)), 24, 24))
%     subplot(1, 4, 2), imagesc(reshape((x_cs_all(:, nn)), 24, 24))
    title('BPDN', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_cs_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_cs_loc(:, nn))), 'FontSize', 18)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)
    
%     subplot(1, 4, 3), imagesc(reshape((x_out_all(:, nn)), 24, 24))
    subplot(1, 4, 3), imagesc(reshape((x_store(:, nn) - x_out_all(:, nn)), 24, 24))
    title('BPDN-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_out_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_out_loc(:, nn))), 'FontSize', 18)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)

%     subplot(1, 4, 4), imagesc(reshape((x_lsm1_all(:, nn)), 24, 24))
    subplot(1, 4, 4), imagesc(reshape((x_store(:, nn) - x_lsm1_all(:, nn)), 24, 24))
    title('RWL1-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_lsm1_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_lsm1_loc(:, nn))), 'FontSize', 18)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)
    drawnow
    set(gcf, 'color', [1,1,1])
    OUT_MOV(nn) = getframe(gcf);
    % pause(0.5)
end
open(OUT_VID);
writeVideo(OUT_VID, OUT_MOV);
close(OUT_VID);

%%

N_thresh = 1000;
thresh_val_cs = logspace(-3, 0, N_thresh);
thresh_val_dcs = logspace(-3, 0, N_thresh);
thresh_val_drw = logspace(-3, 0, N_thresh);

P_MD = zeros(3, N_thresh);
P_FA = zeros(3, N_thresh);

for kk = 1:numel(thresh_val_cs)
    x_loc = x_store ~= 0;
    x_cs_loc = abs(x_cs_all) > thresh_val_cs(kk);
    x_out_loc = abs(x_out_all) > thresh_val_dcs(kk);
    x_lsm1_loc = abs(x_lsm1_all) > thresh_val_drw(kk);

    
    P_MD(:, kk) = [sum(sum(x_loc&(~x_lsm1_loc), 1)), ...
            sum(sum(x_loc&(~x_cs_loc), 1)), ...
        sum(sum(x_loc&(~x_out_loc), 1))];
    
    P_FA(:, kk) = [sum(sum((~x_loc)&x_lsm1_loc, 1)), ...
            sum(sum((~x_loc)&x_cs_loc, 1)), ...
        sum(sum((~x_loc)&x_out_loc, 1))];
end

figure(3)
plot(P_MD.', P_FA.', 'LineWidth', 3)
xlabel('P_{MD}', 'FontSize', 23)
ylabel('P_{FA}', 'FontSize', 23)
set(gca, 'FontSize', 19)


P_tot = P_MD + P_FA;

min_drw = thresh_val_drw(P_tot(1, :) == min(P_tot(1, :)))
min_dcs = thresh_val_dcs(P_tot(3, :) == min(P_tot(3, :)))
min_cs = thresh_val_cs(P_tot(2, :) == min(P_tot(2, :)))

min(P_tot, [], 2)

