%% Mixing layer reconstruction from downsampling
%
% Jared Callaha, 2018

fprintf('Initializing... ')
addpath(genpath('utils'))
load_ml;  % Load and partition mixing layer
nx = flow.nx; ny = flow.ny;
[XX, YY] = meshgrid(x_cc, y_cc);

xratio = 10;  % Sample every so many pixels (total downsampling will be ratio^2)
yratio = 20;


%% Build dictionaries from downsampled training images
Theta = zeros(nx_small*ny_small, mTrain);
for i=1:flow.mTrain
    y = downsample(reshape(Train(:, i), ny, nx), xratio, yratio);
    Theta(:,i) = y(:);
end


%% Build POD library
[U, Sigma, ~] = svd(Train, 'econ'); % Compute POD modes and singular values

% Optimal rank truncation by Gavish-Donoho threshold
beta = size(Train,2)/size(Train,1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(beta,0) * median(diag(Sigma));
r_GD = find(diag(Sigma)>thresh, 1, 'last');

ThetaPCA = zeros(nx_small*ny_small, r_GD);
for i=1:r_GD
    y = downsample(reshape(U(:, i), ny, nx), xratio, yratio);
    ThetaPCA(:, i) = y(:);
end

%% Reconstructions
x = double(Test(:, end));
y = downsample(reshape(x, ny, nx), xratio, yratio);
[ny_small, nx_small] = size(y);
y = y(:);

% Sparse approximation in training library
s_rsr = sp_approx(y, Theta, 0, flow);

% Least-squares estimation of POD coefficients (gappy POD)
cvx_begin
    variable s_pod(r_GD);
    minimize(norm(s_pod, 2));
    subject to
        norm( y-ThetaPCA*s_pod, 2) <= 0.1;
cvx_end;

% Reconstruct fields with both methods
[x_hat, res] = reconstruct(x, Train, full(s_rsr), flow);
[x_pod, res_pod] = reconstruct(x, U(:, 1:r_GD), s_pod, flow);


%% Save reconstructions
save('../output/mixing-layer/downsampled_reconstructions.mat', 'XX', 'YY', 'nx', 'ny', 'flow'...
    , 'x', 'y', 'ny_small', 'nx_small', 'x_hat', 'x_pod', 'res', 'res_pod', 'xratio', 'yratio')


%% Compare POD coefficients of reconstructions
alpha = U(:, 1:r_GD)'*x;
alpha_hat = U(:, 1:r_GD)'*x_hat;
alpha_pod = s_pod;

save('../output/mixing-layer/pod_coeffs.mat', 'alpha', 'alpha_hat', 'alpha_pod')


%% ===========================================
function x_small = downsample(x, xratio, yratio)
    % Window and downsample the 2D snapshot x (for mixing layer)
    [ny, ~] = size(x);

    temp = x(floor(ny/3):floor(2*ny/3) , :);
    x_small = temp(1:yratio:end, 1:xratio:end);
end