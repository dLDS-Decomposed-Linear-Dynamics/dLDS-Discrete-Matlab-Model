
thresh_val_cs = 0.3472; %0.3285;
thresh_val_dcs = 0.2357; %0.2440;
thresh_val_drw = 0.16; %0.385;

x_loc = x_store ~= 0;
x_cs_loc = abs(x_cs_all) > thresh_val_cs;
x_out_loc = abs(x_out_all) > thresh_val_dcs;
x_lsm1_loc = abs(x_lsm1_all) > thresh_val_drw;


cs_image_total = zeros(24, 24, 3);
dcs_image_total = zeros(24, 24, 3);
drw_image_total = zeros(24, 24, 3);

OUT_VID = VideoWriter('Track_Test_Locate.avi');
OUT_VID.FrameRate = 5; %,'Uncompressed AVI');
for nn = 1:size(x_cs_all, 2)
    figure(1)
    subplot(2, 4, 1), imagesc(reshape(x_loc(:, nn), 24, 24))
    title('True Signal', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    cs_image = zeros(24, 24, 3);
    cs_image(:, :, 2) = reshape((~x_cs_loc(:, nn))&x_loc(:, nn), 24, 24);
    cs_image(:, :, 1) = reshape(x_cs_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    cs_image = cs_image + repmat(reshape(x_cs_loc(:, nn)&x_loc(:, nn), 24, 24), [1, 1, 3]);
    cs_image_total = cs_image_total|cs_image;
    
    subplot(2, 4, 2), imagesc(cs_image) %reshape(x_cs_loc(:, nn), 24, 24))
    title('BPDN', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_cs_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_cs_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    
    dcs_image = zeros(24, 24, 3);
    dcs_image(:, :, 2) = reshape((~x_out_loc(:, nn))&x_loc(:, nn), 24, 24);
    dcs_image(:, :, 1) = reshape(x_out_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    dcs_image = dcs_image + repmat(reshape(x_out_loc(:, nn)&x_loc(:, nn), 24, 24), [1, 1, 3]);
    dcs_image_total = dcs_image_total|dcs_image;
    
    subplot(2, 4, 3), imagesc(dcs_image) %reshape(x_out_loc(:, nn), 24, 24))
    title('BPDN-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_out_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_out_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    
    drw_image = zeros(24, 24, 3);
    drw_image(:, :, 2) = reshape((~x_lsm1_loc(:, nn))&x_loc(:, nn), 24, 24);
    drw_image(:, :, 1) = reshape(x_lsm1_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    drw_image = drw_image + repmat(reshape(x_lsm1_loc(:, nn)&x_loc(:, nn), 24, 24), [1, 1, 3]);
    drw_image_total = drw_image_total|drw_image;
    
    subplot(2, 4, 4), imagesc(drw_image) %reshape(x_lsm1_loc(:, nn), 24, 24))
    title('RWL1-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_lsm1_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_lsm1_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    
    subplot(2, 4, [5,6]), plot(cumsum(sum(x_loc(:, 1:nn)&(~x_lsm1_loc(:, 1:nn)), 1)), '-b', 'LineWidth', 3)
    subplot(2, 4, [5,6]), hold on
    plot(cumsum(sum(x_loc(:, 1:nn)&(~x_out_loc(:, 1:nn)), 1)), '--m', 'LineWidth', 3)
    plot(cumsum(sum(x_loc(:, 1:nn)&(~x_cs_loc(:, 1:nn)), 1)), ':k', 'LineWidth', 3)
    
    set(gca, 'XLim', [1, 100], 'YLim', [0, 575])
    title('Missed Detections', 'FontSize', 25)
    set(gca, 'FontSize', 20)
    hold off
    
    subplot(2, 4, [7,8]), plot(cumsum(sum((~x_loc(:, 1:nn))&x_lsm1_loc(:, 1:nn), 1)), '-b', 'LineWidth', 3)
    subplot(2, 4, [7,8]), hold on
    plot(cumsum(sum((~x_loc(:, 1:nn))&x_out_loc(:, 1:nn), 1)), '--m', 'LineWidth', 3)
    plot(cumsum(sum((~x_loc(:, 1:nn))&x_cs_loc(:, 1:nn), 1)), ':k', 'LineWidth', 3)
    set(gca, 'XLim', [1, 100], 'YLim', [0,125])
    title('False Alarms', 'FontSize', 25)
    set(gca, 'FontSize', 20)
    legend('RWL1-DF', 'BPDN-DF', 'BPDN')
    legend('Location', 'NorthWest')
    hold off
    
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

color_limits = [min([min(min(min(x_store))), min(min(min(x_store - x_cs_all))), ...
    min(min(min(x_store - x_out_all))), min(min(min(x_store - x_lsm1_all)))]), ...
    max([max(max(max(x_store))), max(max(max(x_store - x_cs_all))), ...
    max(max(max(x_store - x_out_all))), max(max(max(x_store - x_lsm1_all)))])];

OUT_VID = VideoWriter('Track_Test_Error_withRMSE.avi');
OUT_VID.FrameRate = 5; %,'Uncompressed AVI');
for nn = 1:size(x_cs_all, 2)
    figure(1)
    subplot(2, 4, 1), imagesc(reshape(x_store(:, nn), 24, 24))
    title('True Signal', 'FontSize', 25)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)
    
    subplot(2, 4, 2), imagesc(reshape((x_store(:, nn) - x_cs_all(:, nn)), 24, 24))
%     subplot(1, 4, 2), imagesc(reshape((x_cs_all(:, nn)), 24, 24))
    title('BPDN', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_cs_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_cs_loc(:, nn))), 'FontSize', 18)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)
    
%     subplot(1, 4, 3), imagesc(reshape((x_out_all(:, nn)), 24, 24))
    subplot(2, 4, 3), imagesc(reshape((x_store(:, nn) - x_out_all(:, nn)), 24, 24))
    title('BPDN-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_out_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_out_loc(:, nn))), 'FontSize', 18)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)

%     subplot(1, 4, 4), imagesc(reshape((x_lsm1_all(:, nn)), 24, 24))
    subplot(2, 4, 4), imagesc(reshape((x_store(:, nn) - x_lsm1_all(:, nn)), 24, 24))
    title('RWL1-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_lsm1_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_lsm1_loc(:, nn))), 'FontSize', 18)
    colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [], 'CLim', color_limits)
    
    subplot(2, 4, [5,6,7,8]), plot(sum((x_store(:, 1:nn) - x_lsm1_all(:, 1:nn)).^2, 1)./(sum(x_store(:, 1:nn).^2, 1)), '-b', 'LineWidth', 3)
    subplot(2, 4, [5,6,7,8]), hold on
    plot(sum((x_store(:, 1:nn) - x_out_all(:, 1:nn)).^2, 1)./(sum(x_store(:, 1:nn).^2, 1)), '--m', 'LineWidth', 3)
    plot(sum((x_store(:, 1:nn) - x_cs_all(:, 1:nn)).^2, 1)./(sum(x_store(:, 1:nn).^2, 1)), ':k', 'LineWidth', 3)
    set(gca, 'XLim', [1, 100], 'YLim', [0,0.25])
    title('rMSE', 'FontSize', 25)
    set(gca, 'FontSize', 20)
    legend('RWL1-DF', 'BPDN-DF', 'BPDN')
    legend('Location', 'NorthWest')
    hold off
    
    drawnow
    set(gcf, 'color', [1,1,1])
    OUT_MOV(nn) = getframe(gcf);
    % pause(0.5)
    
end
open(OUT_VID);
writeVideo(OUT_VID, OUT_MOV);
close(OUT_VID);

%%

thresh_val_cs = 0.24;
thresh_val_dcs = 0.139;
thresh_val_drw = 0.72;

x_loc = x_store ~= 0;
x_cs_loc = abs(x_cs_all) > thresh_val_cs;
x_out_loc = abs(x_out_all) > thresh_val_dcs;
x_lsm1_loc = abs(x_lsm1_all) > thresh_val_drw;



cs_image_total = zeros(24, 24, 3);
cs_image_fa = zeros(24, 24, 1);
cs_image_md = zeros(24, 24, 1);
cs_image_cd = zeros(24, 24, 1);
dcs_image_total = zeros(24, 24, 3);
dcs_image_fa = zeros(24, 24, 1);
dcs_image_md = zeros(24, 24, 1);
dcs_image_cd = zeros(24, 24, 1);
drw_image_total = zeros(24, 24, 3);
drw_image_fa = zeros(24, 24, 1);
drw_image_md = zeros(24, 24, 1);
drw_image_cd = zeros(24, 24, 1);

OUT_VID = VideoWriter('Track_Test_Locate_single.avi');
OUT_VID.FrameRate = 5; %,'Uncompressed AVI');
for nn = 1:size(x_cs_all, 2)
    figure(1)
    image_total = image_total|reshape(x_loc(:, nn), 24, 24);
    subplot(1, 4, 1), imagesc(repmat(image_total, [1,1,3]))
    title('True Signal', 'FontSize', 25)
    colormap gray
    axis off
    axis image
    cs_image = zeros(24, 24, 3);
    cs_image_fa = cs_image_fa|reshape(x_cs_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    cs_image_md = cs_image_md|reshape((~x_cs_loc(:, nn))&x_loc(:, nn), 24, 24);
    cs_image_cd = cs_image_cd|reshape(x_cs_loc(:, nn)&x_loc(:, nn), 24, 24);
    
    cs_image_total = repmat(reshape(x_cs_loc(:, nn)&x_loc(:, nn), 24, 24), [1, 1, 3]);
    
    cs_image_total(cs_image_fa==1, 3) = 0;
    cs_image_total(cs_image_fa==1, 2) = 0;
    cs_image_total(cs_image_fa==1, 1) = 1;
    cs_image_total(cs_image_md==1, 2) = 1;
    
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
