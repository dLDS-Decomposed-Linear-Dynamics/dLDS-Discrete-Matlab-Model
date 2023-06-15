%%

load Results/chunk_4x_results.mat
rMSE_TMP1 = rMSE_ALL;

load Results/chunk_4x_results2.mat
rMSE_TMP2 = rMSE_ALL;

load Results/chunk_4x_results3.mat
rMSE_TMP3 = rMSE_ALL;

load Results/chunk_4x_results4.mat
rMSE_TMP4 = rMSE_ALL;

load Results/chunk_4x_results5.mat
rMSE_TMP5 = rMSE_ALL;

load Results/chunk_4x_results6.mat
rMSE_TMP6 = rMSE_ALL;

load Results/chunk_WL1_results.mat
rMSE_TMPW = rMSE_ALL;

rMSE_ALL = zeros(size(rMSE_TMP1,1), size(rMSE_TMP1, 2)+1, 6*size(rMSE_TMP1, 3));
rMSE_ALL(:, 1:end-1, 1:4) = rMSE_TMP1;
rMSE_ALL(:, 1:end-1, 5:8) = rMSE_TMP2;
rMSE_ALL(:, 1:end-1, 9:12) = rMSE_TMP3;
rMSE_ALL(:, 1:end-1, 13:16) = rMSE_TMP4;
rMSE_ALL(:, 1:end-1, 17:20) = rMSE_TMP5;
rMSE_ALL(:, 1:end-1, 21:24) = rMSE_TMP6;
rMSE_ALL(:, end, :) = rMSE_TMPW;

%%

vid_avgs = mean(rMSE_ALL, 1);

%%

calc_percent = 1;
calc_mean_av = 1;
bars_to_plot = [1:3,5,9,7];

thresh_val = .01;
means_vals = mean(rMSE_ALL, 1);
high_ind = find(means_vals(:,1,:)>thresh_val);
low_ind = find(means_vals(:,1,:)<=thresh_val);

fprintf('Num high = %d, Num low = %d\n', numel(high_ind), numel(low_ind))

mean_h_pdiff = zeros(1, size(rMSE_TMP1, 2));
mean_l_pdiff = zeros(1, size(rMSE_TMP1, 2));
mean_a_pdiff = zeros(1, size(rMSE_TMP1, 2));

std_h_pdiff = zeros(1, size(rMSE_TMP1, 2));
std_l_pdiff = zeros(1, size(rMSE_TMP1, 2));
std_a_pdiff = zeros(1, size(rMSE_TMP1, 2));

median_h_pdiff = zeros(1, size(rMSE_TMP1, 2));
median_l_pdiff = zeros(1, size(rMSE_TMP1, 2));
median_a_pdiff = zeros(1, size(rMSE_TMP1, 2));

p25_h_pdiff = zeros(1, size(rMSE_TMP1, 2));
p25_l_pdiff = zeros(1, size(rMSE_TMP1, 2));
p25_a_pdiff = zeros(1, size(rMSE_TMP1, 2));

p75_h_pdiff = zeros(1, size(rMSE_TMP1, 2));
p75_l_pdiff = zeros(1, size(rMSE_TMP1, 2));
p75_a_pdiff = zeros(1, size(rMSE_TMP1, 2));


for kk = 1:size(rMSE_ALL, 2)
    if (calc_percent == 1)&&(calc_mean_av == 0)
        diff_highs = reshape((rMSE_ALL(:, kk, high_ind) - rMSE_ALL(:, 4, high_ind))./rMSE_ALL(:, kk, high_ind), [], 1);
        diff_lows = reshape((rMSE_ALL(:, kk, low_ind) - rMSE_ALL(:, 4, low_ind))./rMSE_ALL(:, kk, low_ind), [], 1);
        diff_all = reshape((rMSE_ALL(:, kk, :) - rMSE_ALL(:, 4, :))./rMSE_ALL(:, kk, :), [], 1);
    elseif (calc_percent == 0)&&(calc_mean_av == 0)
        diff_highs = reshape((rMSE_ALL(:, kk, high_ind) - rMSE_ALL(:, 4, high_ind)), [], 1);
        diff_lows = reshape((rMSE_ALL(:, kk, low_ind) - rMSE_ALL(:, 4, low_ind)), [], 1);
        diff_all = reshape((rMSE_ALL(:, kk, :) - rMSE_ALL(:, 4, :)), [], 1);
    elseif (calc_percent == 1)&&(calc_mean_av == 1)
        diff_highs = reshape((vid_avgs(:, kk, high_ind) - vid_avgs(:, 4, high_ind))./vid_avgs(:, kk, high_ind), [], 1);
        diff_lows = reshape((vid_avgs(:, kk, low_ind) - vid_avgs(:, 4, low_ind))./vid_avgs(:, kk, low_ind), [], 1);
        diff_all = reshape((vid_avgs(:, kk, :) - vid_avgs(:, 4, :))./vid_avgs(:, kk, :), [], 1);
    elseif (calc_percent == 0)&&(calc_mean_av == 1)
        diff_highs = reshape((vid_avgs(:, kk, high_ind) - vid_avgs(:, 4, high_ind)), [], 1);
        diff_lows = reshape((vid_avgs(:, kk, low_ind) - vid_avgs(:, 4, low_ind)), [], 1);
        diff_all = reshape((vid_avgs(:, kk, :) - vid_avgs(:, 4, :)), [], 1);
    elseif (calc_percent == 1)&&(calc_mean_av == 2)
        diff_highs = reshape((rMSE_ALL(:, kk, high_ind) - rMSE_ALL(:, 4, high_ind))./repmat(vid_avgs(:, kk, high_ind), [size(rMSE_ALL, 1), 1, 1]), [], 1);
        diff_lows = reshape((rMSE_ALL(:, kk, low_ind) - rMSE_ALL(:, 4, low_ind))./repmat(vid_avgs(:, kk, low_ind), [size(rMSE_ALL, 1), 1, 1]), [], 1);
        diff_all = reshape((rMSE_ALL(:, kk, :) - rMSE_ALL(:, 4, :))./repmat(vid_avgs(:, kk, :), [size(rMSE_ALL, 1), 1, 1]), [], 1);
    elseif (calc_percent == 0)&&(calc_mean_av == 2)
        diff_highs = reshape((rMSE_ALL(:, kk, high_ind) - rMSE_ALL(:, 4, high_ind)), [], 1);
        diff_lows = reshape((rMSE_ALL(:, kk, low_ind) - rMSE_ALL(:, 4, low_ind)), [], 1);
        diff_all = reshape((rMSE_ALL(:, kk, :) - rMSE_ALL(:, 4, :)), [], 1);
    else
        error('Bad parameters!')
    end
    % Means
    mean_h_pdiff(kk) = mean(diff_highs);
    mean_l_pdiff(kk) = mean(diff_lows);
    mean_a_pdiff(kk) = mean(diff_all);
    
    %stds
    std_h_pdiff(kk) = std(diff_highs)/sqrt(numel(diff_highs));
    std_l_pdiff(kk) = std(diff_lows)/sqrt(numel(diff_lows));
    std_a_pdiff(kk) = std(diff_all)/sqrt(numel(diff_all));
    
    % Medians
    median_h_pdiff(kk) = median(diff_highs);
    median_l_pdiff(kk) = median(diff_lows);
    median_a_pdiff(kk) = median(diff_all);
    
    % 25% and 75%
    high_sort = sort(diff_highs);
    low_sort = sort(diff_lows);
    all_sort = sort(diff_all);
    
    p25_h_pdiff(kk) = high_sort(round(0.25*numel(high_sort)));
    p25_l_pdiff(kk) = low_sort(round(0.25*numel(low_sort)));
    p25_a_pdiff(kk) = all_sort(round(0.25*numel(all_sort)));
    
    p75_h_pdiff(kk) = high_sort(round(0.75*numel(high_sort)));
    p75_l_pdiff(kk) = low_sort(round(0.75*numel(low_sort)));
    p75_a_pdiff(kk) = all_sort(round(0.75*numel(all_sort)));
    
end


mean_alls = (1 + 99*calc_percent)*[mean_a_pdiff(bars_to_plot); mean_h_pdiff(bars_to_plot);mean_l_pdiff(bars_to_plot)];
median_alls = (1 + 99*calc_percent)*[median_a_pdiff(bars_to_plot); median_h_pdiff(bars_to_plot);median_l_pdiff(bars_to_plot)];

pud_all = zeros(3, numel(bars_to_plot), 2);
pud_all(:, :, 1) = median_alls - (1 + 99*calc_percent)*[p25_a_pdiff(bars_to_plot); p25_h_pdiff(bars_to_plot); p25_l_pdiff(bars_to_plot)];
pud_all(:, :, 2) = (1 + 99*calc_percent)*[p75_a_pdiff(bars_to_plot); p75_h_pdiff(bars_to_plot); p75_l_pdiff(bars_to_plot)] - median_alls;


std_alls = zeros(3, numel(bars_to_plot), 2);
std_alls(:, :, 1) = (1 + 99*calc_percent)*[std_a_pdiff(bars_to_plot); std_h_pdiff(bars_to_plot); std_l_pdiff(bars_to_plot)];
std_alls(:, :, 2) = (1 + 99*calc_percent)*[std_a_pdiff(bars_to_plot); std_h_pdiff(bars_to_plot); std_l_pdiff(bars_to_plot)];


%

% % subplot(1, 2, 1), bar(100*[mean_a_pdiff([1:3,5:7]); mean_h_pdiff([1:3,5:7]);mean_l_pdiff([1:3,5:7])]);
% subplot(1, 2, 1), barwitherr(std_alls, mean_alls);
% set(gca, 'FontSize', 20, 'XTick', [1,2,3], 'XTickLabel', {'Mean', 'High Error', 'Low Error'})
% ylabel('Mean % Improvement', 'FontSize', 24)
% % legend('BPDN', 'RWL1', 'BPDN-DF', 'DCS-AMP', 'RWL1-mult', 'modCS')
% % my_xticklabels([1,2,3], {{'Mean %'; 'Improvement'} {'Mean High Error';'Improvement'} {'Mean Low Error';'Improvement'}}, 'FontSize', 20);
% 
% % subplot(1, 2, 2), bar(100*[median_a_pdiff([1:3,5:7]); median_h_pdiff([1:3,5:7]);median_l_pdiff([1:3,5:7])])
% subplot(1, 2, 2), barwitherr(pud_all, median_alls)
% set(gca, 'FontSize', 20, 'XTick', [1,2,3], 'XTickLabel', {'Median', 'High Error', 'Low Error'})
% ylabel('Median % Improvement', 'FontSize', 24)
% % legend('BPDN', 'RWL1', 'BPDN-DF', 'DCS-AMP', 'RWL1-mult', 'modCS')
% % my_xticklabels([1,2,3], {{'Median %'; 'Improvement'} {'Median High Error';'Improvement'} {'Median Low Error';'Improvement'}}, 'FontSize', 20);

%

subplot(1, 2, 1), err_bars_cross = zeros(2, numel(bars_to_plot), 2);
subplot(1, 2, 1), hold on
err_bars_cross(1, :, :) = std_alls(1, :, :);
err_bars_cross(2, :, :) = pud_all(1, :, :);

barwitherr(err_bars_cross, 100*[mean_a_pdiff(bars_to_plot);median_a_pdiff(bars_to_plot)])
set(gca, 'FontSize', 20, 'XTick', [1,2], 'XTickLabel', {'Mean', 'Median'}, 'YLim', [-15,100], 'View', [90, 90])
ylabel('% Improvement', 'FontSize', 24)
subplot(1, 2, 1), hold off

subplot(1, 2, 2), err_bars_cross = zeros(2, numel(bars_to_plot), 2);
subplot(1, 2, 2), hold on
err_bars_cross(1, :, :) = std_alls(2, :, :);
err_bars_cross(2, :, :) = pud_all(2, :, :);

barwitherr(err_bars_cross, 100*[mean_h_pdiff(bars_to_plot);median_h_pdiff(bars_to_plot)])
set(gca, 'FontSize', 20, 'XTick', [1,2], 'XTickLabel', {'Mean', 'Median'}, 'View', [90, 90])
ylabel('% Improvement', 'FontSize', 24)
legend('BPDN', 'RWL1', 'BPDN-DF', 'DCS-AMP', 'WL1P', 'modCS')
subplot(1, 2, 1), hold off

