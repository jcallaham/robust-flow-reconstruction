%% Plot centerline reconstruction examples for the mixing layer and sparsity histogram
%
% Jared Callaham, 2018

close all;
% Load output from mixing-layer/ml_centerline.m
load('output/mixing-layer/centerline_reconstructions.mat')


%% Plot reconstructions and dotted boxes around windows
figure(1)
pcolor(reshape(x, ny, nx)); hold on; axis off
title('Test field')

figure(2)
pcolor(reshape(x+noise, ny, nx)); hold on
ylim([-1 ny+2]); axis off
title('Noisy test field')

figure(3)
pcolor(reshape(x_sr, ny, nx)); hold on; axis off
title(sprintf('Sparse representation (global)\n err = %0.2f', res_sr));

figure(4);
pcolor(reshape(x_pod, ny, [])); hold on; axis off
title(sprintf('Least-squares POD (global)\n err = %0.2f', res_pod))

figure(5)
pcolor(reshape(x_kernel_sr, ny, nx)); hold on; axis off
title(sprintf('Sparse representation (kernels), k=%d\nerr=%0.2f', k, res_kernel_sr));

figure(6)
pcolor(reshape(x_kernel_pod, ny, nx)); hold on; axis off
title(sprintf('Least squares POD (kernels), k=%d\nerr=%0.2f', k, res_kernel_pod));

figure(7);
pcolor(reshape(x_window_sr, ny, [])); hold on; axis off
title('Sparse representation (window)')

figure(8);
pcolor(reshape(x_window_pod, ny, [])); hold on; axis off
title('Least-squares POD (window)')


% Plot dashed line around windows
for i = 1:k
    for j=[2 7 8]
        figure(j)
        plot([x0(i) x0(i)], [0 ny], 'k--', 'LineWidth', 1.5)
        plot([x0(i) x0(i+1)], [ny, ny], 'k--', 'LineWidth', 1.5)
        plot([x0(i+1) x0(i)], [0 0], 'k--', 'LineWidth', 1.5)
        if i==length(x0)-1
            plot([x0(i+1) x0(i+1)], [ny 0], 'k--', 'LineWidth', 1.5)
        end
    end
end

% Plot sensor locations
sensor_locs = [floor(sensor_idx/ny); mod(sensor_idx, ny)];
for i=2:8
    figure(i);
    scatter(sensor_locs(1, :), sensor_locs(2, :), 20, 'k', 'Filled')
end

% Format and save
labels = {'a', 'b', 'c', 'd', 'e', 'f'};
for i=1:8
    figure(i)
    colormap(flow.cmap); shading flat; caxis(flow.clim); axis on
    set(gcf, 'Position', [100*i 200*i 960, 225])
    
    if SAVE_FIGS && i<7  % Last two figures didn't go in the paper
        name = sprintf('figs/FIG8_%s.png', labels{i});
        saveas(gcf, name);
    end
end

%% Plot histogram of sparsity for local vs global reconstruction
load('output/mixing-layer/sparsity.mat')

figure();
bin_edges = 0:1e-4:0.05;

histogram(K_global_sr, 20, 'Normalization', 'probability'); hold on
histogram(K_window_sr, 10, 'Normalization', 'probability');

xlabel('Relative sparsity')
ylabel('Frequency')

grid on
set(gcf, 'Position', [100 100 800 200])

if SAVE_FIGS
    saveas('FIG8_hist.svg', gcf)
end