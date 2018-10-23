close all;

fprintf('Initializing... \n')
addpath(genpath('../utils'))
load_ml;  % Load and partition mixing layer
nx = flow.nx; ny = flow.ny;
[XX, YY] = meshgrid(x_cc, y_cc);

xratio = 10;  % Sample every so many pixels (total downsampling will be ratio^2)
yratio = 10;

energy_rescale = true;


%% Construct linearly growing windows as a function of downstream coordinate
window_width = @(x) floor(0.5*(x/nx)*ny) + 20;
x0 = 1;
window_end = x0 + window_width(x0);  % Keep track of where the last window is
while window_end < nx  % Note that last one will be too far (stop the reconstruction loop short)
    x0 = [x0 window_end]; 
    window_end = window_end + window_width(window_end);
end
% Extend last window to end of domain
x0(end) = nx;
k = length(x0)-1;  % Number of kernels/windows


%% Define kernels
kernel_x = floor(x0(1:end-1) + 0.5*diff(x0));  % Center points of kernels - same as windows
kernel_y = ny/2 + 0*kernel_x; % Center on midline

kernel_width = nx/(2*k);  % Width of kernels
phi = @(x) exp(-x.^2/kernel_width^2);

%% Compute normalized kernels
denom = zeros(ny, nx);
for j=1:k
    [XX, YY] = meshgrid((1:nx) - kernel_x(j), (1:ny) - kernel_y(j) );
    RR = sqrt(XX.^2 + YY.^2);  % Compute distance from every point to kernel
    Phi{j} = phi(RR);
    denom = denom + Phi{j};
end

% Normalize and reshape to diagonal
for j=1:k
    kernel = Phi{j}./denom;
    kernel(kernel < 1e-2) = 0;
    Phi{j} = spdiags(kernel(:), 0, flow.n, flow.n);
end
fprintf('Done\n')

%% Compute downsampling parameters
row_idx = floor(2*ny/5):yratio:floor(3*ny/5); % Middle 1/5 of flow=
col_idx = 1:xratio:nx;

% Compute sampled points
[CIDX, RIDX] = meshgrid(col_idx, row_idx);
sensor_idx = CIDX(:)*ny + RIDX(:);  % Convert to 1D indices
[ny_small, nx_small] = size(RIDX);

%% Plot example downsampling
x = double(Test(:, end));
y = x(sensor_idx);
sensor_locs = [floor(sensor_idx/ny) mod(sensor_idx, ny)];


%% Build libraries
% Library from downsampled training images
D = double(Train(sensor_idx, :));

% POD library
[U, Sigma, ~] = svd(Train, 'econ');
beta = size(Train,2)/size(Train,1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(beta,0) * median(diag(Sigma));
r_GD = find(diag(Sigma)>thresh, 1, 'last');
r = 50;

%% Reconstruct
disp('Global reconstructions...')
x = double(Test(:, end));

eta = 0;
noise = eta*flow.avg_energy*randn(size(x));

y = x(sensor_idx)+noise(sensor_idx);

s_rsr = sp_approx(y, D, eta, flow);
s_pod = U(sensor_idx, 1:r)\y;

[x_hat, res] = reconstruct(x, Train, full(s_rsr), flow, energy_rescale);
[x_pod, res_pod] = reconstruct(x, U(:, 1:r), s_pod, flow, energy_rescale);

%% Window reconstruction
disp('Window reconstructions...')

x_window_sr = zeros(size(x));
x_window_pod = zeros(size(x));
x_kernel = zeros(size(x));

for j =1:k   % Loop through kernels/windows
   
    %% Instantaneous window reconstruction
    window =  ny*x0(j):ny*(x0(j+1))-1;  % Window for comparison
    window_sensors = sensor_idx(ismember(sensor_idx, window)); % Find sensors in window
    D = double(Train(window_sensors, :));
    y = x(window_sensors) + noise(window_sensors);
    
    local_flow = flow;
    local_flow.avg_energy = mean( sqrt( mean(Train(window, :).^2, 1) ) );
    
    %% Sparse approximation in training library
    s = sp_approx(y, D, eta, flow);  
    [window_hat, ~] = reconstruct(x(window), Train(window, :), s, local_flow, energy_rescale);
    
    x_window_sr(window) = window_hat;
    
    
    %% L2 POD reconstruction in windows
    % Local POD modes - use fewer than num sensors
    r_window = floor(0.5*length(y));
    [U_window, ~, ~] = svd(Train(window, :), 'econ');
    D = zeros(flow.n, r_window);
    D(window, :) = U_window(:, 1:r_window);
    s = D(window_sensors, :)\y;
    [window_hat, ~] = reconstruct(x(window), U_window(:, 1:r_window), s, local_flow, energy_rescale);
    x_window_pod(window) = window_hat;
    
    %% Kernel sparse reconstruction
    % "kernelized" measurement matrix, library, and state (to save memory)
    kernel_idx = find(diag(Phi{j}));
    kernel = Phi{j}(kernel_idx, kernel_idx);
    kernel_sensors = sensor_idx(ismember(sensor_idx, kernel_idx)); % Restrict to sensors in kernel
    
    C = full(Phi{j}(kernel_sensors, :));  % measurement matrix
    library = kernel*double(Train(kernel_idx, :));
    x_j = kernel*x(kernel_idx);
    
    y = C*(x + noise);
    D = double(C*Train);
    
    local_flow.avg_energy = mean( sqrt( mean(library.^2, 1) ) );
    
    s = sp_approx(y, D, eta, flow);
    [kernel_hat, ~] = reconstruct(x_j, library, s, local_flow, true);
    x_kernel(kernel_idx) = x_kernel(kernel_idx) + kernel_hat;

end

res_window_sr = norm(x_window_sr - x)/norm(x);  % Sparse approximation error
res_window_pod = norm(x_window_pod - x)/norm(x);  % L2 POD error
res_kernel = norm(x_kernel - x)/norm(x);
    

%% Save reconstructions
save('../output/mixing-layer/downsampled_reconstructions.mat', 'eta', 'XX', 'YY', 'nx', 'ny', 'flow'...
    , 'x', 'y', 'ny_small', 'nx_small', 'x_hat', 'x_pod', 'res', 'res_pod', 'xratio', 'yratio'...
    , 'x_window_sr', 'x_window_pod', 'x_kernel')
