% Robust reconstruction of sea surface temperature fields
% NOAA OISST v2 data set
%
% Jared Callaham 2018

warning('off', 'all');
addpath(genpath('../utils'))

%% Load and partition data
% Note: no mean subtraction (already done)
load_sst;
nx = flow.nx; ny = flow.ny;

energy_rescale = false;

%% Compute POD
[U, Sigma, ~] = svd(Train, 'econ');  % Singular value decomposition
beta = size(Train, 2)/size(Train, 1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(beta,0) * median(diag(Sigma));
r_GD = find(diag(Sigma)>thresh, 1, 'last');

%% Reconstruction from random points
t = 835;  % Choose test field (something far from training set)
% Choose randomly from valid points between 50 deg S and 50 deg N
ns = 10;
sensor_idx = choose_sensors(ns, flow.mask, ny, lat);

rand_locs = [mod(sensor_idx, ny) round(sensor_idx/ny)];  % (x, y) coordinates for plotting

% Point measurement
x = Test(:, t);              % Choose test field
y = double(x(sensor_idx)); % Noisy slice

% Sparse approximation using training dictionary
D = double(Train(sensor_idx, :));
s = var_approx(y, D, 0.5);  % Optimal relaxation
% Reconstruct temperature field (without energy rescaling)
[x_sr, res_sr] = reconstruct(x, Train, full(s), flow, energy_rescale); 
disp(res_sr);

%% Least-squares POD with same measurements
r = 2;
alpha = U(sensor_idx, 1:r)\y;
[x_pod, res_pod] = reconstruct(x, U(:, 1:r), alpha, flow, false);
disp(res_pod)

%% Reconstruction with Gappy POD (using many more points)
% Choose random sensor locoations
ns = 500;
sensor_idx = choose_sensors(ns, flow.mask, ny, lat);
gappy_locs = [mod(sensor_idx, ny) round(sensor_idx/ny)];  % (x, y) coordinates for plotting

% Optimal sensor placement (via QR-conditioning)
%[~, ~, sensor_idx] = qr(Ur(:, 1:min(ns, r_GD))', 'vector');
%sensor_idx = sensor_idx(1:ns);

% Gappy POD reconstruction
y = x(sensor_idx);
r = 50;
alpha = U(sensor_idx, 1:r)\y;
[x_gappy, res_gappy] = reconstruct(x, U(:, 1:r), alpha, flow, false);
disp(res_gappy);

%% Save results for plotting
x(flow.mask) = NaN;
x_sr(flow.mask) = NaN;
x_pod(flow.mask) = NaN;
x_gappy(flow.mask) = NaN;

save('../output/noaa-sst/reconstruction_examples.mat', 'x', 'x_sr', 'x_pod'...
    , 'x_gappy', 'rand_locs', 'gappy_locs', 'res_sr', 'res_pod', 'res_gappy', 'flow')

%% Choose random, but appropriate, sensors
function sensor_idx = choose_sensors( ns, mask, ny, lat )
% Choose random valid measurement locations on (50S, 50N)
n = length(mask);
sensor_idx = zeros(ns, 1);
for i=1:ns
    found = false;
    while(~found)
        sensor_idx(i) = randi(n);
        yloc = lat(mod(sensor_idx(i), ny)+1);  % Latitude corresponding to y-coordinat in grid
        found = ~mask(sensor_idx(i)) && (yloc > -50) && (yloc < 50);
    end
end

end

% Same as sp_approx function, but with variable relaxation epsilon
function s = var_approx(y, D, eps)
    m = size(D, 2);
    cvx_begin quiet;
        variable s(m);
        minimize( norm(s,1) );
        subject to
            norm(D*s - y, 2) <= eps;
    cvx_end;

    s = double(s);  % Convert from cvx variable to regular float
end

