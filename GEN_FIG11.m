%% Plot example HYCOM reconstructions
% Jared Callaham 2018

close all;
clearvars -except SAVE_FIGS
load('data/gom_vort.mat', 'lat', 'lon')
% Load output from hycom/hycom_reconstruct.m
load('output/hycom/reconstruction_examples.mat')

%% Plot results
mask = flow.mask;
mean_flow(mask) = NaN;  % Mask NaN values (so we can color land differently)
x_measure(mask) = NaN;
x_measure(~mask & x_measure~=0) = 1;  % Just plot black dots for measurements

% NaN mask for plotting background
x(mask) = NaN;
x_sr(mask) = NaN;
x_pod(mask) = NaN;
x_local_sr(mask) = NaN;
x_local_pod(mask) = NaN;

ny = length(lat); nx = length(lon);


%%
figure()
% Plot test field
pcolor(reshape(x+flow.mean_flow, ny, nx))
shading interp; caxis(flow.clim);
title('Test field')
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])
colormap(flow.cmap)

% Plot measured points
figure()  % Save axis to tweak colormap later
pcolor(reshape(x_measure, ny, nx))
shading flat; caxis([0, 1]);
title('Measurement locations')
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])
colormap(1-gray);

% Reconstructed field from sparse representation 
figure()
pcolor(reshape(x_sr+flow.mean_flow, ny, nx))
shading interp; caxis(flow.clim);
title(sprintf('Sparse reconstruction\n ns = 4k, err: %0.2f', res_sr))
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])
colormap(flow.cmap)

% Reconstructed field by gappy POD
figure()
pcolor(reshape(x_pod+flow.mean_flow, ny, nx))
shading interp; caxis(flow.clim);
title(sprintf('SVHT Gappy POD (r=2760, ns=4k)\n err: %0.2f', res_pod))
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])
colormap(flow.cmap)

% Reconstruction from kernel sparse representation
figure()
pcolor(reshape(x_local_sr+flow.mean_flow, ny, nx))
shading interp; caxis(flow.clim);
title(sprintf('Kernel sparse reconstruction (k=%d, ns=4k)\n err: %0.2f', k, res_local_sr))
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])
colormap(flow.cmap)

% Reconstruction from kernel gappy POD
figure()
pcolor(reshape(x_local_pod+flow.mean_flow, ny, nx))
shading interp; caxis(flow.clim);
title(sprintf('Kernel POD reconstruction (k=%d, ns=4k)\n err: %0.2f', k, res_local_pod))
set(gca, 'Color', [100 100 100]/255);  % Set NaN (as background) to gray
set(gca,'xtick',[]); set(gca,'ytick',[])
colormap(flow.cmap)

%%
labels = {'a', 'b', 'c', 'd', 'e', 'f'};
for i=1:6
    set(figure(i) ,'Position',[100*i 100*i 2*length(lon) 2*length(lat)])
    
    % Note: may have to save manually to get background
    if SAVE_FIGS
        saveas(gcf, sprintf('figs/FIG11%s.png', labels{i}))
    end
end


