%% Plot SST anomaly reconstructions
%
% Jared Callaham

clearvars -except SAVE_FIGS
addpath(genpath('utils'));

% Load output from noaa-sst/sst_reconstruct.m
run('noaa-sst/sst_reconstruct.m')
load('output/noaa-sst/reconstruction_examples.mat')

%% Plot commands
figure()

% Plot truth
subplot(221)
pcolor(reshape(x, flow.ny, flow.nx));
shading interp, caxis(flow.clim)
title('Snapshot')
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])

% Plot sparse reconstruction in training library
subplot(222)
pcolor(reshape(x_sr, flow.ny, flow.nx)); hold on
shading interp, caxis(flow.clim)
scatter(rand_locs(:, 2), rand_locs(:, 1), 25, 'k', 'filled')
title(sprintf('%d measurements\n err: %0.4f', size(rand_locs, 1), res_sr) )
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])

% Plot least-squares POD
subplot(223)
pcolor(reshape(x_pod, flow.ny, flow.nx)); hold on
shading interp, caxis(flow.clim)
scatter(rand_locs(:, 2), rand_locs(:, 1), 25, 'k', 'filled')
title(sprintf('Gappy POD\n err: %0.4f', res_pod) )
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])

% Plot gappy POD (many more measurements)
subplot(224)
pcolor(reshape(x_gappy, flow.ny, flow.nx)); hold on
shading interp, caxis(flow.clim)
scatter(gappy_locs(:, 2), gappy_locs(:, 1), 25, 'k', 'filled')
title(sprintf('%d measurements (gappy POD)\n err: %0.4f', size(gappy_locs, 1), res_gappy) )
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])

colormap(flow.cmap);
set(gcf,'Position',[100 100 6*flow.nx 6*flow.ny])

%% Note: may have to save manually to get background (land) to be gray
if SAVE_FIGS
    saveas(gcf, 'figs/FIG9.svg')
end
