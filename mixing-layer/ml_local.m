%% Mixing layer reconstruction with downstream windowing

fprintf('Initializing... ')
addpath(genpath('utils'))
load_ml;  % Load and partition mixing layer
nx = flow.nx; ny = flow.ny;
[XX, YY] = meshgrid(x_cc, y_cc);

%% Construct linearly growing windows as a function of downstream coordinate

% Automatically build linearly growing windows
window_width = @(x) floor(0.5*(x/nx)*ny);  % Linearly growing window width

x0 = round(0.14*nx);
window_end = x0 + window_width(x0);  % Keep track of where the last window is
while window_end < nx  % Note that last one will be too far (stop the reconstruction loop short)
    x0 = [x0 window_end]; 
    window_end = window_end + window_width(window_end);
end
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

beta = 0.5*size(Train,2)/size(Train,1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(beta,0) * median(diag(Sigma));
r_GD = find(diag(Sigma)>thresh, 1, 'last');
Ur = double(U(:, 1:r_GD));  % Truncate
clearvars U % Free up memory
fprintf('Done\n')

%% Reconstruction
disp('Reconstructing...')
t = 100; % Choose test snapshot
x = double(Test(:, t));

sigma = 0.7;  % Noise level
rms_vort = double(mean( sqrt( mean(Train.^2, 1) ) ));  % For noise level only - no energy rescaling
noise = sigma*rms_vort*randn([flow.n, 1]);

x_window_sr = zeros(size(x));
x_window_pod = zeros(size(x));
x_kernel_sr = zeros(size(x));
x_kernel_pod = zeros(size(x));

% Find window and global sensor locations
sensor_idx = []; % All sensor locations (for later global reconstruction)
for i=1:k
    % Stagger sensors along midline of window
    x_span = floor( (0:ns-1)*window_width(x0(i))/ns );  % Evenly space sensors along midline of window
    sensor_locs = [x_span + x0(i); ny/2+0*x_span]';   % Horizontal slice
    sensor_idx = [sensor_idx (sensor_locs(:, 1)*ny + sensor_locs(:, 2))'];  % Save sensor indices for global reconstruction
end

for j =1:k   % Loop through kernels/windows
   
    %% Instantaneous window reconstruction
    window =  ny*x0(j):ny*(x0(j+1))-1;  % Window for comparison
    window_sensors = sensor_idx(ismember(sensor_idx, window)); % Find sensors in window
    D = double(Train(window_sensors, :));
    y = x(window_sensors) + noise(window_sensors);
    eps = 3*sigma*rms_vort*sqrt(length(y));  % Tolerance for optimization problem
   
    s = sp_approx(y, D, sigma, flow);
    % Reconstruct (zeros for mean, since this wasn't subtracted), NaN for
    % energy rescaling (not justified in this case?)
    [window_hat, res] = reconstruct(x(window), Train(window, :), s, flow);
    
    x_window_sr(window) = window_hat;
    
    
    %% L2 POD reconstruction in windows
    s = lsqminnorm(Ur(window_sensors, :), y);
    [window_hat, pod_res] = reconstruct(x(window), Ur(window, :), s, flow);
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
    
    s = sp_approx(y, D, sigma, flow);
    [window_hat, ~] = reconstruct(x_j, library, s, flow);
    x_kernel_sr(kernel_idx) = x_kernel_sr(kernel_idx) + window_hat;
    
    %% Kernel POD reconstruction
    s = lsqminnorm(C*Ur, y);
    [window_hat, ~] = reconstruct(x_j, kernel*Ur(kernel_idx, :), s, flow);
    x_kernel_pod(kernel_idx) = x_kernel_pod(kernel_idx) + window_hat;
    
    fprintf('%d/%d  res:%0.2f   POD:%0.2f\n', j, k, res, pod_res);
end

% Compute total residuals for kernels
res_kernel_sr = norm(x_kernel_sr - x)/norm(x);
res_kernel_pod = norm(x_kernel_pod - x)/norm(x);

%% Global reconstruction
D = double(Train(sensor_idx, :));
y = x(sensor_idx) + noise(sensor_idx);

s = sp_approx(y, D, sigma, flow);
% Reconstruct (zeros for mean, since this wasn't subtracted)
[x_sr, res_sr] = reconstruct(x, Train, s, flow);
disp(res)

%% L2 POD
s = lsqminnorm(Ur(sensor_idx, :), y);
[x_pod, res_pod] = reconstruct(x, Ur, s, flow);
disp(res_pod)



%% Save results

save('../output/mixing-layer/local_reconstructions.mat', 'flow', 'x', 'noise', 'nx', 'ny'...
    , 'k', 'x0', 'x_window_sr', 'x_window_pod', 'x_kernel_sr', 'x_kernel_pod', 'x_sr'...
    , 'x_pod', 'sensor_idx', 'res_sr', 'res_pod', 'res_kernel_sr', 'res_kernel_pod')
    
