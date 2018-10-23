%% Reconstruction of HYCOM data ('92-'18)
% Relatively sparse measurement (~10%): compare sparse reconstruction with
% OMP to gappy POD (least squares)
%
% Jared Callaham 2018

clear all
addpath(genpath('../utils'))

%% Load vorticity data
disp('Initializing...')
load_hycom;
nx = flow.nx; ny = flow.ny; mask = flow.mask;
n = nx*ny;

%% Calculate POD of training data
[U, Sigma, ~] = svd(double(Train), 'econ');
Sigma = diag(Sigma);
% Gavish-Donoho optimal hard threshold
beta = size(Train,2)/size(Train,1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(beta,0) * median(Sigma);
r_GD = find(Sigma>thresh, 1, 'last');
Ur = U(:, 1:r_GD);
clearvars U


%% Construct dictionaries
disp('Building dictionaries')

ns = 4000;  % Number of sensors
sensor_idx = randsample(find(~mask), ns);  % Choose randomly from valid locations
t = 50;
x = double(Test(:, t));  
y = x(sensor_idx);  % Essentially C*x

%%
% Form a full-state vector but only nonzero entries are measurement (just for visualization)
x_measure = zeros(n, 1);
x_measure(sensor_idx) = y+flow.mean_flow(sensor_idx);

% Dictionary from measurements of training set 
D = double(Train(sensor_idx, :));  
    
%% Reconstruction - compare test image, global sparse reconstruction, gappy POD, and local sparse reconstruction

%% OMP
% omp_rescale = zeros(size(Train, 2), 1);
% % Normalize dictionary columns
% for j=1:flow.mTrain
%     omp_rescale(j) = norm(D(:, j));
%     D(:, j) = D(:, j)/omp_rescale(j);
% end

%K = 400;    % desired sparsity of solution (empirically estimated)
% % Identify sparse coefficients with matching pursuit
%s_sr = omp(D, y, D'*D, K); 
%[x_sr, res_sr] = reconstruct(x, Train, full(s_sr).\omp_rescale, flow);

%% L1  (much slower than OMP)
disp('Computing sparse approximation')
% Identify sparse coefficients
s_sr = sp_approx(y, D, 0, flow);
[x_sr, res_sr] = reconstruct(x, Train, full(s_sr), flow);
% Reconstruct from sparse representation
K_sr = sum(abs(s_sr)>1e-6);

disp('Computing gappy POD approximation')
% Least-squares POD (global)
s_pod = (double(Ur(sensor_idx, :)))\y;  % Coefficients from least squares solution
        
% Reconstruction
[x_pod, res_pod] = reconstruct(x, Ur, s_pod, flow);
K_pod = sum(abs(s_pod) > 1e-6);

%% Local sparse reconstruction
disp('Computing local reconstructions')

% Define and compute normalized kernels
kx = 12; ky = 8;  % Number of kernels to use
k = kx*ky;

kernel_x = floor(linspace(1, nx, kx));  % Center points of kernels
kernel_y = floor(linspace(1, ny, ky));  % Center points of kernels
[kernel_x, kernel_y] = meshgrid(kernel_x, kernel_y);
kernel_x = kernel_x(:); kernel_y = kernel_y(:);


sigma = nx/(2*kx);  % Width of kernels
phi = @(x) exp(-x.^2/sigma^2);  % kernel function
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
    Phi{j} = spdiags(kernel(:), 0, n, n);
end

% Solve each reconstruction problem independently
x_local_sr = zeros(n, 1);
x_local_pod = zeros(n, 1);
K_local_sr = [];
K_local_pod = [];
for j=1:k
    % "kernelized" measurement matrix, library, and state (to save memory)
    kernel_idx = find(diag(Phi{j}));
    kernel = Phi{j}(kernel_idx, kernel_idx);
    kernel_sensors = sensor_idx(ismember(sensor_idx, kernel_idx)); % Restrict to sensors in kernel
    
    C = full(Phi{j}(kernel_sensors, :));  % measurement matrix
    library = kernel*double(Train(kernel_idx, :));  % Reweighted library for reconstruction
    x_j = kernel*x(kernel_idx);   % Reweighted state in this kernel (for error comparison)
    local_flow = flow;
    local_flow.mean_flow = kernel*flow.mean_flow(kernel_idx);
    
    y = C*x;
    if ~isempty(y)  % If kernel isn't just land
        fprintf('%d/%d: %d sensors\n', j, k, length(y));
        
        % Sparse representation
        D = double(C*Train);
        s = sp_approx(y, D, 0, flow);
        [x_hat, ~] = reconstruct(x_j, library, s, local_flow);
        x_local_sr(kernel_idx) = x_local_sr(kernel_idx) + x_hat;
        K_local_sr = [K_local_sr sum(abs(s)>1e-6)];  % Compute sparsity of solution
        
        % L2 POD
        s = (C*Ur)\y;
        [x_hat, ~] = reconstruct(x_j, kernel*Ur(kernel_idx, :), s, local_flow);
        x_local_pod(kernel_idx) = x_local_pod(kernel_idx) + x_hat;
        K_local_pod = [K_local_pod sum(abs(s)>1e-6)];  % Compute sparsity of solution
    end
end

res_local_sr = norm(x_local_sr(~mask) - x(~mask))/norm(x(~mask) + flow.mean_flow(~mask));
res_local_pod = norm(x_local_pod(~mask) - x(~mask))/norm(x(~mask) + flow.mean_flow(~mask));

%%

save('../output/hycom/reconstruction_examples.mat', 'x', 'sensor_idx', 'x_measure' ...
    , 's_sr', 'x_sr', 'res_sr', 'x_pod', 's_pod', 'res_pod', 'flow'...
    , 'x_local_sr', 'res_local_sr', 'x_local_pod', 'res_local_pod'...
    , 'K_sr', 'K_pod', 'K_local_sr', 'K_local_pod', 'k')
    
