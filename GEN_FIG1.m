%% Plot subfigures for diagram of SRC method with cylinder
% Jared Callaham 2018

% Output is partitioned and mean-subtracted Train/Test data and flow params struct
run('re100-cylinder/load_cylinder.m');
nx = flow.nx; ny = flow.ny;
n = nx*ny;

%% Plot snapshots for "dictionary" cartoon
for t=1:7:57
    disp(t);
    
    x = VORTALL(:, t);
    plotCylinder(reshape(x, nx, ny), flow);
    axis off
    if (SAVE_FIGS)
        saveas(gcf, sprintf('figs/FIG1_atom%d.svg', t))
    end
end

%% Plot example figure with point measurements
% Measurement: random pixels
% ns = 5;  % Number of point measurements

% % Restrict to cylinder wake: final 80% of width, middle 50% of height
% sensor_idx = [randperm(round(0.8*ny), ns)' randperm(round(0.5*nx), ns)']; % Choose sensors on restricted area
% sensor_idx = [round(0.2*ny)+sensor_idx(:, 1) round(0.25*nx)+sensor_idx(:, 2)];  % Translate to wake
%  
% % Convert to measurement matrix
% C = spdiags(ones(n, 1), 0, n, n);
% C = C(sensor_idx(:, 2) + (nx-1)*sensor_idx(:, 1), :);  % Sparse measurement matrix

load('output/re100-cylinder/sensor_loc5.mat')   % Predetermined random sensor locations

% Example flow snapshot
x = VORTALL(:, 20);  
plotCylinder(reshape(x, nx, ny), flow); hold on
scatter(sensor_idx(:, 1), sensor_idx(:, 2), 50, 'k', 'filled'); axis off
set(gcf,'Position',[100 100 600 260])
if (SAVE_FIGS)
    saveas(gcf, 'figs/FIG1_sensors.svg')
end

% Show training dictionary
mTrain = 32;
D = C*VORTALL(:, 1:2:mTrain);    % Measured library (includes mean for visualization)
figure()
pcolor([D; D(end, :)]);
colormap(flow.cmap); caxis(flow.clim);
set(gcf,'Position',[1000 100 650 200]); axis off;
if (SAVE_FIGS)
    saveas(gcf, 'figs/FIG1_dict.svg')
end


%% Test image examples
Train = Train(:, 1:2:end);  % Some snapshots are skipped so aspect ratio of library is good
rms_vort = flow.avg_energy;

x = Test(:, 12);
D = C*Train;  % Measured library (now mean-subtracted)

% Add noise to the measurement
eta = 0.65;
noise = eta*rms_vort*randn([n, 1]);
y = C*(x+noise);

% Compute sparse approximation to the entire flow field and calculate error
s = sp_approx(y, D, eta, flow);
[x_hat, res] = reconstruct(x, Train, s, flow);
disp(res)

%% Plot results
figure()
plotCylinder(reshape(x+noise+flow.mean_flow, nx, ny), flow); hold on
scatter(sensor_idx(:, 1), sensor_idx(:, 2), 50, 'k', 'filled'); axis off
set(gcf,'Position',[100 100 600 260]);
title('y')
if (SAVE_FIGS)
    saveas(gcf, 'figs/FIG1_x.svg')
end

figure()
pcolor([y y; y(end) y(end)]);
colormap(flow.cmap); caxis(flow.clim);
set(gcf,'Position',[1000 100 39 175]); axis off;
title('y')
if (SAVE_FIGS)
    saveas(gcf, 'figs/FIG1_y.svg')
end

s_plot = 1-s;  % Invert for plotting
figure()
pcolor([s_plot s_plot; s_plot(end) s_plot(end)]);
set(gcf,'Position',[1000 100 38 730]); axis off;
title('s')
colormap gray;
if (SAVE_FIGS)
    saveas(gcf, 'figs/FIG1_s.svg')
end

figure()
plotCylinder(reshape(x_hat+flow.mean_flow, nx, ny), flow); axis off
set(gcf,'Position',[100 100 600 260]);
title('xhat')

if (SAVE_FIGS)
    saveas(gcf, 'figs/FIG1_xhat.svg')
end
