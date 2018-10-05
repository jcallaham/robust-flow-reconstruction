%% Plot examples of different measurement strategies along with err vs noise/corruption

clearvars -except SAVE_FIGS
addpath(genpath('utils'))
addpath('re100-cylinder')
run('re100-cylinder/load_cylinder.m');
run('re100-cylinder/define_measurements.m');  % Generate or load measurement matrices

nx = flow.nx; ny = flow.ny;

%% Plot greyscale snapshot
x = Test(:, 17);   % Test field
vort_field = abs(x + flow.mean_flow); 
vort_field = 1 - (vort_field/max(vort_field));  % Normalize to (0, 1)
plotCylinder(reshape(vort_field, nx, ny), flow);
colormap(gray)
caxis([.83, 1]);
load('utils/colorscheme.mat')  % From colorbrewer2.org - for plot lines

%% Plot measurement types
% Window
y_lo = min(window_loc(:, 1));  y_hi = max(window_loc(:, 1));
x_lo = min(window_loc(:, 2));  x_hi = max(window_loc(:, 2));
p = fill([y_lo y_hi y_hi y_lo], [x_lo x_lo x_hi x_hi], colorscheme(1, :));
p.FaceAlpha = 0.3;
p.EdgeColor = colorscheme(1, :);
p.LineWidth = 3;

% Horizontal slice
plot(hslice_loc(:, 1), hslice_loc(:, 2), 'Color', colorscheme(3, :), 'LineWidth', 3)

% Vertical slice
plot(vslice_loc(:, 1), vslice_loc(:, 2), 'Color', colorscheme(4, :), 'LineWidth', 3)

% Random points
scatter(rand_loc(:, 1), rand_loc(:, 2), 70, colorscheme(2, :), 'filled');

axis off;
if SAVE_FIGS
	saveas(gcf, 'figs/FIG5a.svg')
end


%% Plot reconstruction error vs corruption (figure 5c)
load('output/re100-cylinder/corrupt_loop_out.mat') % Output from re100-cylinder/corrupt_loop.m
figure(); hold on;
for i=1:4
    errorbar(rho, res(i, :), err(i, :), 'o', 'Color', colorscheme(i, :)...
        ,'MarkerSize', 5, 'MarkerFaceColor', colorscheme(i, :), 'LineStyle', '--'...
        ,'LineWidth', 1.5)
end
xlabel('Corruption percentage \rho')
ylabel('Reconstruction error')

grid on; box on;

set(gcf,'Position',[100 100 600 260])
if SAVE_FIGS
	saveas(gcf, 'figs/FIG5c.svg')
end


%% Plot reconstruction error vs noise (figure 5e)
load('output/re100-cylinder/noise_loop_out.mat')  % Output of re100-cylinder/noise_loop.m
figure(); hold on;
for i=1:4
    errorbar(sigma, res(i, :), err(i, :), 'o', 'Color', colorscheme(i, :)...
        ,'MarkerSize', 5, 'MarkerFaceColor', colorscheme(i, :), 'LineStyle', '--')
end
xlabel('Noise level \sigma')
ylabel('Reconstruction error')
ylim([0 .2])
grid on;

box on;
set(gcf,'Position',[100 100 600 260])
if SAVE_FIGS
	saveas(gcf, 'figs/FIG5e.svg')
end


%% Plot test field (figure 5b)
figure()
plotCylinder(reshape(x, nx, ny), flow); axis off
title('Test field')

if SAVE_FIGS
	saveas(gcf, 'figs/FIG5b.svg')
end


%% Plot example reconstruction from vertical slice (figure 5d/f)
sigma = 0.2;  % Dense noise level
rho = 0;  % Percentage of corrupt grid points
% Reconstruct cylinder flow field in separate script
[x_noisy, x_hat, res] = cylinder_reconstruct(x, C_vslice, Train, flow, sigma, rho);    
figure()
subplot(121)
plotCylinder(reshape(x_noisy+flow.mean_flow,nx,ny), flow); hold on, axis off
colorbar();

% Vertical slice
plot(vslice_loc(:, 1), vslice_loc(:, 2), 'Color', 'k', 'LineWidth', 3)
title('Corrupted snapshot');

subplot(122)
plotCylinder(reshape(x_hat+flow.mean_flow, nx, ny), flow);
axis off
title(sprintf('Training dictionary, L1\n err = %0.4f', res));
set(gcf,'Position',[100 100 1200 260])

if SAVE_FIGS
	saveas(gcf, 'figs/FIG5df.svg')
end
