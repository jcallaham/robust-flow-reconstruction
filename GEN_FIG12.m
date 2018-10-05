%% Plot residual of POD subspace projection versus length of training data
% Jared Callaham 2018

load('output/subspace_residuals.mat')

%% Plot residuals vs m
figure();
subplot(121)
l = plot(m_cyl, cyl_res, '.'); hold on
l.MarkerFaceColor = l.Color;
l = plot(m_mix, mix_res, '.');
l.MarkerFaceColor = l.Color;
l = plot(m_sst, sst_res, '.');
l.MarkerFaceColor = l.Color;
l = plot(m_gom, gom_res, '.');
l = plot(m_gom, rand_gom, '--', 'Color', l.Color);
%plot(m_gom, rand_res);

set(gca, 'XScale', 'Log'); set(gca, 'YScale', 'Log')
xlabel('Training set length m')
ylabel('Error in subspace projection')
legend('Cylinder', 'Mixing Layer', 'SST anomaly', 'HYCOM GoM', 'Random GoM', 'Location', 'Southwest')
ylim([1e-4, 1.5])
grid on;

subplot(122)
plot(m_mix, window_res, 'k.', 'MarkerSize', 4);
grid on; hold on
l = plot(m_mix, mix_res, '.', 'Color', [0.8500    0.3250    0.0980], 'MarkerSize', 4);
set(gca, 'XScale', 'Linear'); set(gca, 'YScale', 'Log');
ylim([1e-4, 1.5]);
xlabel('Training set length m')


set(gcf,'Position',[300 300 590 296])
if SAVE_FIGS
	saveas(gcf, 'figs/FIG12.svg')
end