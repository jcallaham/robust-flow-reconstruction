%% Plot reconstruction error vs number of sensors for all flows

close all;
clearvars -except SAVE_FIGS
load('output/err_vs_ns.mat') % Load output of err_vs_sensors.m

figure()

%% Plot errorbar lines
l_cyl = errorbar(ns, cyl_res(1, :), cyl_err(1, :), 'Marker', 'o'); hold on;
l_mix = errorbar(ns, ml_res(1, :), ml_err(1, :), 'Marker', 'o'); hold on;
l_sst = errorbar(ns, sst_res, sst_err, 'Marker', 'o'); hold on;
l_gom = errorbar(ns, gom_res, gom_err, 'Marker', 'o'); hold on;

%% Cylinder with noise
l_cyl.MarkerFaceColor = l_cyl.Color;
l_cyl.LineWidth = 2;
l_cyl.MarkerSize = 5;
l_cyl.LineStyle = '--';

% plot(ns, cyl_res(2, :), 'Color',  l_cyl.Color, 'LineWidth', 1.5);
% plot(ns, cyl_res(3, :), 'Color',  l_cyl.Color, 'LineWidth', 1);
% plot(ns, cyl_res(4, :), 'Color',  l_cyl.Color, 'LineWidth', 0.5);

%% Mixing layer with noise
l_mix.MarkerFaceColor = l_mix.Color;
l_mix.LineWidth = 2;
l_mix.MarkerSize = 5;
l_mix.LineStyle = '--';

% plot(ns, ml_res(2, :), 'Color',  l_mix.Color, 'LineWidth', 1.5);
% plot(ns, ml_res(3, :), 'Color',  l_mix.Color, 'LineWidth', 1);
% plot(ns, ml_res(4, :), 'Color',  l_mix.Color, 'LineWidth', 0.5);

%% SST
l_sst.MarkerFaceColor = l_sst.Color;
l_sst.LineWidth = 2;
l_sst.MarkerSize = 5;
l_sst.LineStyle = '--';

%% HYCOM
l_gom.MarkerFaceColor = l_gom.Color;
l_gom.LineWidth = 2;
l_gom.MarkerSize = 5;
l_gom.LineStyle = '--';

%%
set(gca, 'XScale', 'log')
grid on;
ylim([0 2.5])

legend('Cylinder', 'Mixing Layer', 'SST', 'GoM')
xlabel('# random point measurements')
ylabel('Reconstruction error')

set(gcf,'Position',[100 100 800 400])
if SAVE_FIGS
	saveas(gcf, 'figs/FIG13.svg')
end