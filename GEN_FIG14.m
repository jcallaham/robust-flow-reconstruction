close all;
load('output/hycom/kernel_viz.mat')

%% Plot weighted sum of approximate ranks
nx = flow.nx; ny = flow.ny; mask = flow.mask;
n = nx*ny;

% Plot un-normalized kernel coverage
unweighted_kernels(mask) = NaN;

mean_flow = flow.mean_flow;
mean_flow(mask) = NaN; % Mask land for plotting

figure()
subplot(121)
pcolor(reshape(mean_flow, ny, nx)); hold on
colormap(gca, flow.cmap); caxis(flow.clim); colorbar; shading interp

% Plot contours for unweighted kernels
contour(reshape(unweighted_kernels, ny, nx), [0.7, 0.8, 0.9, 0.95], 'Color', 'k');
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])
title('Mean vorticity with kernel centers')


% POD 99% energy cutoff
weighted_dim(mask) = NaN;

subplot(122)
pcolor(reshape(weighted_dim, ny, nx)); shading flat
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
colormap(gca, 'jet'); colorbar
set(gca,'xtick',[]); set(gca,'ytick',[])
title('Local "rank" (99% energy cutoff)')

set(gcf ,'Position',[100 100 4*nx 2*ny])
