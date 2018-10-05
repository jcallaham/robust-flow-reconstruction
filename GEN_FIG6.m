%% Plot reconstruction error vs noise for different reconstruction strategies
% Jared Callaham 2018

load('output/re100-cylinder/reconst_loop_out.mat') % Output from re100-cylinder/reconstruction_loop.m
load('utils/colorscheme.mat')  % From colorbrewer2.org - for plot lines

%% Make errorbar plot
figure(); hold on; grid on;
for i=1:size(res, 1)
    l = errorbar(eta, res(i, :), err(i, :));
    l.Color = colorscheme(i, :);
    l.LineStyle = '--';
    l.LineWidth=1.5;
    l.Marker = markers(i);
    l.MarkerSize=10;
    l.MarkerFaceColor = colorscheme(i, :);
end
    
legend('Training, L1', 'K-SVD, L1', 'POD, L1', 'POD, L2', 'Training, L2', 'Location', 'NorthWest')
xlim([0, 0.2])
xlabel('Noise level')
ylabel('Reconstruction residual')

box on
set(gcf,'Position',[100 100 800 400])
if SAVE_FIGS
	saveas(gcf, 'figs/FIG6.svg')
end

%% Generate insets as examples
% Load and partition cylinder data
run('re100-cylinder/load_cylinder.m')
mTrain = size(Train, 2);

load('output/re100-cylinder/measurement_matrices.mat', 'C_vslice');

rms_vort = flow.avg_energy;
C = C_vslice;


x = Test(:, 18);  % Choose test field
sigma = 0.08;  % Noise level (scaled by rms vorticity)
noise = sigma*rms_vort*randn(size(x));
eps = 3*sigma*rms_vort*sqrt(size(C, 1));  % Tolerance for this noise level 

y = C*(x + noise);  % Noisy measurement
D = double(C*Train);  % Training library
% Solve L2 optimization problem in training library
cvx_begin;
    variable s(mTrain);
    minimize(norm(s, 2));
    subject to
        norm(D*s - y, 2) <= eps;
cvx_end;
[x_hat, r] = reconstruct(x, Train, s, flow);

figure();
subplot(311)
plotCylinder(reshape(x_hat+flow.mean_flow, nx, ny), flow); axis off;
title(sprintf('Training dictionary, L2\n sigma = %0.2f, err = %0.4f', sigma, r));

% Reconstruction by sparse approximation in the training library
s = sp_approx(y, D, sigma, flow);
[x_hat, r] = reconstruct(x, Train, s, flow);

subplot(312)
plotCylinder(reshape(x_hat+flow.mean_flow, nx, ny), flow)
axis off
title(sprintf('Training dictionary, L1\n sigma = %0.2f, err = %0.4f', sigma, r));

sigma = 0.2; % Noise level (fraction of rms vorticity)
noise = sigma*rms_vort*randn(size(x));
eps = 3*sigma*rms_vort*sqrt(size(C, 1));  % Tolerance for this noise level 

y = C*(x + noise); % Noisy measurement

[U, Sigma, ~] = svd(Train, 'econ');  % Compute library of POD modes

% Optimal rank truncation of POD library (by Gavish-Donoho threshold)
bsigma = size(Train, 2)/size(Train, 1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(bsigma,0) * median(diag(Sigma));
m_pod = find(diag(Sigma)>thresh, 1, 'last'); % Rank truncation with Gavish-Donoho SVHT
ThetaPCA = U(:, 1:m_pod)*Sigma(1:m_pod, 1:m_pod);

D_pod = C*ThetaPCA;  % Measured POD library

% Reconstruction by sparse approximation in POD library
s = sp_approx(y, D_pod, sigma, flow);
[x_hat, r] = reconstruct(x, ThetaPCA, s, flow);

subplot(313)
plotCylinder(reshape(x_hat+flow.mean_flow, nx, ny), flow)
axis off
title(sprintf('POD, L1\n sigma = %0.2f, err = %0.4f', sigma, r));


if SAVE_FIGS
	saveas(gcf, 'figs/FIG6_insets.svg')
end