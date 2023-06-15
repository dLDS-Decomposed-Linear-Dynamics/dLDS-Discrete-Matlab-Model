%%

thresh_val_cs = 0.24;
thresh_val_dcs = 0.139;
thresh_val_drw = 0.72;

x_loc = x_store ~= 0;
x_cs_loc = abs(x_cs_all) > thresh_val_cs;
x_out_loc = abs(x_out_all) > thresh_val_dcs;
x_lsm1_loc = abs(x_lsm1_all) > thresh_val_drw;

image_total = zeros(24, 24, 1);
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
    
    cs_image_fa = cs_image_fa|reshape(x_cs_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    cs_image_md = cs_image_md|reshape((~x_cs_loc(:, nn))&x_loc(:, nn), 24, 24);
    cs_image_cd = cs_image_cd|reshape(x_cs_loc(:, nn)&x_loc(:, nn), 24, 24);
    
    cs_image_total = repmat(cs_image_cd, [1, 1, 3]);
    
    TMP = cs_image_total(:, :, 3);
    TMP(cs_image_fa==1) = 0;
    TMP(cs_image_md==1) = 0;
    cs_image_total(:, :, 3) = TMP;
    TMP = cs_image_total(:, :, 2);
    TMP(cs_image_fa==1) = 0;
    TMP(cs_image_md==1) = 1;
    cs_image_total(:, :, 2) = TMP;
    TMP = cs_image_total(:, :, 1);
    TMP(cs_image_md==1) = 0;
    TMP(cs_image_fa==1) = 1;
    cs_image_total(:, :, 1) = TMP;
    
    
    
    subplot(1, 4, 2), imagesc(cs_image_total) %reshape(x_cs_loc(:, nn), 24, 24))
    title('BPDN', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_cs_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_cs_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    
    dcs_image_fa = dcs_image_fa|reshape(x_out_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    dcs_image_md = dcs_image_md|reshape((~x_out_loc(:, nn))&x_loc(:, nn), 24, 24);
    dcs_image_cd = dcs_image_cd|reshape(x_out_loc(:, nn)&x_loc(:, nn), 24, 24);
    
    dcs_image_total = repmat(dcs_image_cd, [1, 1, 3]);
    
    TMP = dcs_image_total(:, :, 3);
    TMP(dcs_image_fa==1) = 0;
    TMP(dcs_image_md==1) = 0;
    dcs_image_total(:, :, 3) = TMP;
    TMP = dcs_image_total(:, :, 2);
    TMP(dcs_image_fa==1) = 0;
    TMP(dcs_image_md==1) = 1;
    dcs_image_total(:, :, 2) = TMP;
    TMP = dcs_image_total(:, :, 1);
    TMP(dcs_image_md==1) = 0;
    TMP(dcs_image_fa==1) = 1;
    dcs_image_total(:, :, 1) = TMP;
    
    
    subplot(1, 4, 3), imagesc(dcs_image_total) %reshape(x_out_loc(:, nn), 24, 24))
    title('BPDN-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_out_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_out_loc(:, nn))), 'FontSize', 18)
%     colormap gray
%     axis off
    axis image
    set(gca, 'XTick', [], 'YTick', [])
    
    drw_image_fa = drw_image_fa|reshape(x_lsm1_loc(:, nn)&(~x_loc(:, nn)), 24, 24);
    drw_image_md = drw_image_md|reshape((~x_lsm1_loc(:, nn))&x_loc(:, nn), 24, 24);
    drw_image_cd = drw_image_cd|reshape(x_lsm1_loc(:, nn)&x_loc(:, nn), 24, 24);
    
    drw_image_total = repmat(drw_image_cd, [1, 1, 3]);
    
    TMP = drw_image_total(:, :, 3);
    TMP(drw_image_fa==1) = 0;
    TMP(drw_image_md==1) = 0;
    dcs_image_total(:, :, 3) = TMP;
    TMP = drw_image_total(:, :, 2);
    TMP(drw_image_fa==1) = 0;
    TMP(drw_image_md==1) = 1;
    drw_image_total(:, :, 2) = TMP;
    TMP = drw_image_total(:, :, 1);
    TMP(drw_image_md==1) = 0;
    TMP(drw_image_fa==1) = 1;
    drw_image_total(:, :, 1) = TMP;
    
    
    subplot(1, 4, 4), imagesc(drw_image_total) %reshape(x_lsm1_loc(:, nn), 24, 24))
    title('RWL1-DF', 'FontSize', 25)
%     xlabel(sprintf('P(md) = %d, P(fa) = %d', sum(x_loc(:, nn)&(~x_lsm1_loc(:, nn))), ...
%         sum((~x_loc(:, nn))&x_lsm1_loc(:, nn))), 'FontSize', 18)
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
