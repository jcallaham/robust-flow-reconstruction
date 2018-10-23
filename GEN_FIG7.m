%% Plot downsampled reconstruction examples for the mixing layer
%
% Jared Callaham, 2018

clearvars -except SAVE_FIGS
close all;
% Load output from mixing-layer/ml_downsample.m
load('output/mixing-layer/downsampled_reconstructions.mat')


%% Plot results
figure(1)
pcolor(XX, YY, reshape(x, ny, nx)); hold on;
colorbar();
title('Test field')

figure(2);
pcolor(reshape(y, ny_small, nx_small))
title('Downsampled')

figure(3)
pcolor(XX, YY, reshape(x_hat, ny, nx));
title(sprintf('Sparse reconstruction: err = %0.2f', res))

figure(4)
pcolor(XX, YY, reshape(x_pod, ny, nx));
title(sprintf('Gappy POD: err = %0.2f', res_pod))


%% Format and save
labels = {'a', 'b', 'c', 'd'};
for i=1:4
    figure(i)
    colormap(flow.cmap); shading flat; caxis(flow.clim);
    set(gcf, 'Position', [100*i 200*i 960, 225])
    
    if i==2
        set(gcf, 'Position', [100*i 200*i 960 75])
    end
    
    if SAVE_FIGS
        name = sprintf('figs/FIG7_%s.png', labels{i});
        saveas(gcf, name);
    end
end


%% Plot histogram of sparsity for local vs global reconstruction
load('output/mixing-layer/downsampled_sparsity.mat')

figure();
bin_edges = 0:1e-4:0.05;

histogram(K_global_sr, 20, 'Normalization', 'probability'); hold on
histogram(K_window_sr, 10, 'Normalization', 'probability');

xlabel('Relative sparsity')
ylabel('Frequency')

grid on
set(gcf, 'Position', [100 100 800 200])

if SAVE_FIGS
    saveas('FIG7_hist.svg', gcf)
end