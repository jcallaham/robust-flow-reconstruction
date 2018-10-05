%% Plot singular value spectrum for all data sets
% Jared Callaham 2018

load('output/sv_spectra.mat')

%% Plot singular value spectra
figure();
subplot(211);
plot(cyl_sv/sum(cyl_sv), '.'); hold on;
plot(mix_sv/sum(mix_sv), '.');
plot(sst_sv/sum(sst_sv), '.');
plot(gom_sv/sum(gom_sv), '.');

set(gca, 'XScale', 'Log'); set(gca, 'YScale', 'Log')
grid on;
xlabel('POD mode')
ylim([1e-6 1])
ylabel('Normalized singular value')
legend('Cylinder', 'Mixing Layer', 'SST anomaly', 'HYCOM GoM')

subplot(212);
plot(cumsum(cyl_sv)/sum(cyl_sv), '.'); hold on;
plot(cumsum(mix_sv)/sum(mix_sv), '.');
plot(cumsum(sst_sv)/sum(sst_sv), '.');
plot(cumsum(gom_sv)/sum(gom_sv), '.');

set(gca, 'XScale', 'Log');
grid on;
xlabel('POD mode')
ylabel('Cumulative energy')

set(gcf,'Position',[300 300 590 650])
if SAVE_FIGS
	saveas(gcf, 'figs/FIG3.svg')
end
