clear all; close all; clc

addpath(genpath('../utils'));
load_hycom;

nx = flow.nx; ny = flow.ny; mask = flow.mask;
n = nx*ny;

Train(mask, :) = 0; % For computing SVD

%% Define kernels
kx = 12; ky = 8; % Number of kernels to use
k = kx*ky;

kernel_x = floor(linspace(1, nx, kx));  % Center points of kernels
kernel_y = floor(linspace(1, ny, ky));  % Center points of kernels
[kernel_x, kernel_y] = meshgrid(kernel_x, kernel_y);
kernel_x = kernel_x(:); kernel_y = kernel_y(:);


sigma = nx/(2*kx);  % Width of kernels
phi = @(x) exp(-x.^2/sigma^2);  % kernel function

% % Get rid of kernels centered on land
% mask_2D = reshape(mask, ny, nx);
% kernel_x = kernel_x( diag(~mask2D(kernel_y, kernel_x))); % Get rid of masked locations

%% Compute normalized kernels
denom = zeros(ny, nx);
for j=1:kx*ky
    [XX, YY] = meshgrid((1:nx) - kernel_x(j), (1:ny) - kernel_y(j) );
    RR = sqrt(XX.^2 + YY.^2);  % Compute distance from every point to kernel
    Phi{j} = phi(RR);
    denom = denom + Phi{j};
end
 


%% Normalize and reshape to diagonal
for j=1:kx*ky
    kernel = Phi{j}./denom;
    kernel(kernel < 1e-2) = 0;
    Phi{j} = spdiags(kernel(:), 0, n, n);
end


%% Break up kernels and calculate emprical "rank"
weighted_dim = zeros(n, 1);
for j=1:k
    % "kernelized" measurement matrix, library, and state (to save memory)
    kernel_idx = find(diag(Phi{j}));
    kernel = Phi{j}(kernel_idx, kernel_idx);
    
    if ~all(mask(kernel_idx))  % Check that there are non-land values
        library = kernel*double(Train(kernel_idx, :));  % Reweighted library for reconstruction
        [~, Sigma, ~] = svd(library, 'econ');
        sv = diag(Sigma)/sum(diag(Sigma));

        r = find( cumsum(sv)>0.99, 1 )-1; % Find number of POD modes to capture 99% energy
        weighted_dim(kernel_idx) = weighted_dim(kernel_idx) + diag(kernel)*r;
        fprintf('%d/%d: r=%d\n', j, k, r);
    end
end

%% Plot weighted local sparsity
load('../output/hycom/reconstruction_examples.mat', 'K_local_sr')
K = K_local_sr;
weighted_K = zeros(n, 1);
idx = 1; % Use to skip empty kernels
for j=1:k
    % "kernelized" measurement matrix, library, and state (to save memory)
    kernel_idx = find(diag(Phi{j}));
    kernel = Phi{j}(kernel_idx, kernel_idx);
    
    if ~all(mask(kernel_idx))  % Check that there are non-land values
        weighted_K(kernel_idx) = weighted_K(kernel_idx) + diag(kernel)*K(idx);
        fprintf('%d/%d: K=%d\n', j, k, K(idx));
        
        idx = idx+1;
    end
end

%%
unweighted_kernels = denom(:);
save('../output/hycom/kernel_viz.mat', 'flow', 'unweighted_kernels'...
    , 'weighted_dim', 'weighted_K')
