%%
load Results/2013_01_25_Psweep_40t_m70.mat
Fin_Mat2_TMP = Fin_Mat2;
load Results/2013_01_26_Psweep_40t_m70_s20_05xVar
Fin_Mat2_TMP4 = Fin_Mat2;
load Results/2013_01_26_Psweep_40t_m70_s25
Fin_Mat2_TMP3 = Fin_Mat2;
load Results/2013_01_26_Psweep_40t_m70_s20_2xVar
Fin_Mat2_TMP2 = Fin_Mat2;
load Results/2013_01_26_Psweep_40t_m70_s15

LL_pt = 1;
UL_pt = 8;

figure(1), hold off
figure(1)
plot(20*poi_vec(LL_pt:UL_pt)*2, Fin_Mat2_TMP(LL_pt:UL_pt, 4), 's--c', 'LineWidth', 3)
figure(1), hold on
plot(20*poi_vec(LL_pt:UL_pt)*2, Fin_Mat2_TMP2(LL_pt:UL_pt, 4), '^:b', 'LineWidth', 3)
plot(20*poi_vec(LL_pt:UL_pt)*2, Fin_Mat2_TMP4(LL_pt:UL_pt, 4), 'o-k', 'LineWidth', 3)
plot(25*poi_vec(LL_pt:UL_pt)*2, Fin_Mat2_TMP3(LL_pt:UL_pt, 4), '.-.r', 'LineWidth', 3)
plot(15*poi_vec(LL_pt:UL_pt)*2, Fin_Mat2(LL_pt:UL_pt, 4), 'v-.g', 'LineWidth', 3)
figure(1), hold off;

    
legend('S=20,\sigma^2 = 0.001', 'S=20,\sigma^2 = 0.002', 'S=20,\sigma^2 = 0.0005', 'S=25,\sigma^2 = 0.001', 'S=15,\sigma^2 = 0.001')
xlabel('Innovations Sparsity', 'FontSize', 30)
ylabel('Steady State rMSE', 'FontSize', 30)
set(gca, 'FontSize', 23, 'Ylim', [0,0.3], 'Xlim', ...
    [s_num*poi_vec(LL_pt)*2, s_num*poi_vec(UL_pt)*2])

%% Histogram plot

rMSE_ALL2 = rMSE_ALL(:, [1:8], :);

FULL_MEANS = mean(mean(rMSE_ALL2, 3), 1);
FULL_MEDIANS = median(mean(rMSE_ALL2, 3), 1);
FIG_HEIGHT = 150;
Y_LIMS_PER1 = [70, 75, 35, 140, 90, 75, 55, 60];
alg_names = {'BPDN', 'RWL1', 'BPDN-DF', 'RWL1-DF', 'DCS-AMP', 'RWL1-SS', 'modCS', 'WL1P'};
bin_cents = linspace(min(min(rMSE_ALL2)), 0.1, 30);
figure(1);
for kk = 1:8
    subplot(2, 4, kk), hold on
    hist(mean(rMSE_ALL2(:, kk, :), 3), bin_cents)
    line([FULL_MEANS(kk), FULL_MEANS(kk)], [0, FIG_HEIGHT], 'LineStyle', ...
        '--', 'LineWidth', 3, 'Color', [0,0,0])
    box on
    title(alg_names{kk}, 'FontSize', 24)
    set(gca, 'FontSize', 20, 'YLim', [0, Y_LIMS_PER1(kk)], ...
        'XLim', [0, 0.065], 'XTick', [0.02,0.04,0.06])
    xlabel('rMSE', 'FontSize', 24)
    axis square
    plotarrow('Start', [FULL_MEDIANS(kk), -15*Y_LIMS_PER1(kk)/150], 'Stop', ...
        [FULL_MEDIANS(kk), 0], 'Length', 4, 'BaseAngle', 85, ...
        'TipAngle', 35, 'Width', 3)
    hold off
end