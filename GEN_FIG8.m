%% Plot estimated POD coefficients for downsampled reconstruction of the mixing layer

load('output/mixing-layer/pod_coeffs.mat')


figure()
plot(alpha, 'k:', 'LineWidth', 2); hold on   % Actual coefficients from projecting test field
plot(alpha_hat, '-', 'LineWidth', 2);  % Coefficients by projecting sparse approximation
plot(alpha_pod, '--', 'LineWidth', 2);  % Coefficients estimated by direct L2 regression

xlim([0, 50])
grid on;
xlabel('POD mode')
ylabel('Coefficient')
legend('Flow field', 'Sparse representation', 'Gappy POD')
ylim([-6.5, 6.5])

set(gcf, 'Position', [300 300 800 300])

%%
if SAVE_FIGS
    saveas(gcf, 'figs/FIG8.svg')
end