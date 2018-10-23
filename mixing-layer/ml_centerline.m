%% Mixing layer reconstruction with downstream windowing

fprintf('Initializing... ')
addpath(genpath('../utils'))
%load_ml;  % Load and partition mixing layer
nx = flow.nx; ny = flow.ny;
[XX, YY] = meshgrid(x_cc, y_cc);

energy_rescale = true;  % Rescale for reconstructions

%% Construct linearly growing windows as a function of downstream coordinate

% Automatically build linearly growing windows
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
ns = 10;  % Evenly spaced along midline


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

%% POD library
fprintf('Computing POD... ')
[U, Sigma, ~] = svd(Train, 'econ');
r = 50;
Ur = double(U(:, 1:r));  % Truncate
clearvars U % Free up memory
fprintf('Done\n')

%% Reconstruction
disp('Reconstructing...')
t = 240; % Choose test snapshot (100 is OK)
x = double(Test(:, t));

eta = 0.3;  % Noise level
noise = eta*flow.avg_energy*randn([flow.n, 1]);

x_window_sr = zeros(size(x));
x_window_pod = zeros(size(x));
x_kernel_sr = zeros(size(x));
x_kernel_pod = zeros(size(x));

% Find window and global sensor locations
sensor_idx = []; % All sensor locations (for later global reconstruction)
for i=1:k
    % Stagger sensors along midline of window
    x_span = floor( (1:ns)*(x0(i+1)-x0(i))/(ns+1) );  % Evenly space sensors along midline of window
    sensor_locs = [x_span + x0(i); ny/2+0*x_span]';   % Horizontal slice
    sensor_idx = [sensor_idx (sensor_locs(:, 1)*ny + sensor_locs(:, 2))'];  % Save sensor indices for global reconstruction
end

r_window = 8;  % Number of local modes to keep
for j =1:k   % Loop through kernels/windows
   
    %% Instantaneous window reconstruction
    window =  ny*x0(j):ny*(x0(j+1))-1;  % Window for comparison
    window_sensors = sensor_idx(ismember(sensor_idx, window)); % Find sensors in window
    D = double(Train(window_sensors, :));
    y = x(window_sensors) + noise(window_sensors);
    local_flow = flow;
    local_flow.avg_energy = mean( sqrt( mean(Train(window, :).^2, 1) ) );
    
    s = sp_approx(y, D, eta, flow);
    % Reconstruct
    [window_hat, res] = reconstruct(x(window), Train(window, :), s, local_flow, energy_rescale);
    
    x_window_sr(window) = window_hat;
    
    %% L2 POD reconstruction in windows
    % Local POD modes - use fewer than num sensors
    [U_window, ~, ~] = svd(Train(window, :), 'econ');
    D = zeros(flow.n, r_window);
    D(window, :) = U_window(:, 1:r_window);
    s = D(window_sensors, :)\y;
    [window_hat, pod_res] = reconstruct(x(window), U_window(:, 1:r_window), s, local_flow, energy_rescale);
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
    [window_hat, ~] = reconstruct(x_j, library, s, local_flow, energy_rescale);
    x_kernel_sr(kernel_idx) = x_kernel_sr(kernel_idx) + window_hat;
    
    %% Kernel POD reconstruction
    s = lsqminnorm(C*Ur, y);
    [window_hat, ~] = reconstruct(x_j, kernel*Ur(kernel_idx, :), s, flow, energy_rescale);
    x_kernel_pod(kernel_idx) = x_kernel_pod(kernel_idx) + window_hat;
    
    fprintf('%d/%d  res:%0.2f   POD:%0.2f\n', j, k, res, pod_res);
end

% Compute total residuals for kernels and windows 
res_window_sr = norm(x_window_sr - x)/norm(x);  % Sparse approximation error
res_window_pod = norm(x_window_pod - x)/norm(x);  % L2 POD error
res_kernel_sr = norm(x_kernel_sr - x)/norm(x);
res_kernel_pod = norm(x_kernel_pod - x)/norm(x);

res_global_windowed = norm(x_sr - x)/norm(x);  % Sparse approximation error

%% Global reconstruction
D = double(Train(sensor_idx, :));
y = x(sensor_idx) + noise(sensor_idx);

s = sp_approx(y, D, eta, flow);
% Reconstruct (zeros for mean, since this wasn't subtracted)
[x_sr, res_sr] = reconstruct(x, Train, s, flow, energy_rescale);
disp(res_sr)

%% L2 POD
s = lsqminnorm(Ur(sensor_idx, :), y);
[x_pod, res_pod] = reconstruct(x, Ur, s, flow, energy_rescale);
disp(res_pod)

%% Save results

save('../output/mixing-layer/centerline_reconstructions.mat', 'flow', 'x', 'noise', 'nx', 'ny'...
    , 'k', 'x0', 'x_window_sr', 'x_window_pod', 'x_kernel_sr', 'x_kernel_pod', 'x_sr'...
    , 'x_pod', 'sensor_idx', 'res_sr','res_pod', 'res_kernel_sr', 'res_kernel_pod')
    
