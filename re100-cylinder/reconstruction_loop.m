%% Robust reconstruction of Re=100 cylinder flow
% Test reconstruction error vs noise for different reconstruction strategies
% NOTE: Training set with L2 norm gives same result as PCA L2
%
% Jared Callaham, 2018

warning('off', 'all');
addpath(genpath('../utils'))

load_cylinder;  % Load and partition full data set and dimensions
define_measurements;  % Define or load measurement matrices C

% Partition into mean-subtracted training and test data
rms_vort = flow.avg_energy;
mTrain = size(Train, 2);

C = C_vslice;   % Select from loaded matrices
ns = size(C, 1);

energy_rescale = true;

%% POD library
[U, Sigma, ~] = svd(Train, 'econ');  % Compute modes and singular values

% Optimal rank truncation with Gavish-Donoho threshold
sigma = size(Train, 2)/size(Train, 1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(sigma,0) * median(diag(Sigma));
r_GD = find(diag(Sigma)>thresh, 1, 'last'); % Rank truncation with Gavish-Donoho SVHT
ThetaPCA = U(:, 1:r_GD)*Sigma(1:r_GD, 1:r_GD);

%% K-SVD
% (parameters are hand tuned)
K = 2;
nD = mTrain;

params.data = Train;
params.Tdata = K;
params.dictsize = nD;
params.iternum = 50;
params.memusage = 'high';
[ThetaKSVD,g,~] = ksvd(params,'tr');
save('../output/re100-cylinder/ksvd_dict.mat', 'ThetaKSVD', 'nD')

load('../output/re100-cylinder/ksvd_dict.mat')

%% Loop through noise levels
sigma = 0:0.01:0.2;            % Noise level as a percentage of mean variance
res = zeros(5, length(sigma));  % Normalized residual error
res2 = zeros(5, length(sigma));  % Same, but squared (to calculate variance/error)
num_configs = 1;   % Number of independent noise realizations to test
mTest = size(Test, 2);

% Predefine dictionaries for speed
D_train = C*Train;
D_ksvd = C*ThetaKSVD;
D_pod = C*ThetaPCA;

disp('Beginning loop')
for i=1:length(sigma)  % 1:length(sigma)
    eps = max(3*sigma(i)*rms_vort*sqrt(ns), 0.01);   % Note that L2 methods are sensitive to the tolerance
    for j=1:num_configs
        for t=1:mTest
            x = Test(:, t);
            
            %% Corrupt with noise
            noise = sigma(i)*rms_vort*randn(size(x));
            y = C*(x + noise); % Noisy slice
            
            %% Reconstructions
            % Training dictionary, L1 norm
            s = sp_approx(y, D_train, sigma(i), flow);
            [~, r] = reconstruct(x, Train, s, flow, energy_rescale);
            res(1, i) = res(1, i) + r;
            res2(1, i) = res2(1, i) + r^2;

            % K-SVD dictionary, L1
            s = sp_approx(y, D_ksvd, sigma(i), flow);
            [~, r] = reconstruct(x, ThetaKSVD, s, flow, energy_rescale);
            res(2, i) = res(2, i) + r;
            res2(2, i) = res2(2, i) + r^2;

            % PCA dictionary, L1
            s = sp_approx(y, D_pod, sigma(i), flow);
            [~, r] = reconstruct(x, ThetaPCA, s, flow, energy_rescale);
            res(3, i) = res(3, i) + r;
            res2(3, i) = res2(3, i) + r^2;

            % POD library, L2
            cvx_begin quiet;
                variable s(r_GD);  
                minimize( norm(s,2) );
                subject to
                    norm(D_pod*s - y,2) <= eps;
            cvx_end;
            [~, r] = reconstruct(x, ThetaPCA, s, flow, energy_rescale);
            res(4, i) = res(4, i) + r;
            res2(4, i) = res2(4, i) + r^2;
            
            % Training dictionary, L2
            cvx_begin quiet;
                variable s(mTrain);
                minimize(norm(s, 2));
                subject to
                    norm(D_train*s - y, 2) <= eps;
            cvx_end;
            [~, r] = reconstruct(x, Train, s, flow, energy_rescale);
            res(5, i) = res(5, i) + r;
            res2(5, i) = res2(5, i) + r^2;
        end
    end
    res(:, i) = res(:, i)/(num_configs*mTest);    % Avg error
    res2(:, i) = res2(:, i)/(num_configs*mTest);  % Avg error, squared
    fprintf('Noise level: \t%0.2f\n', sigma(i))
    fprintf('Training (L1): \t%0.4f\n', res(1, i));
    fprintf('K-SVD (L1): \t%0.4f\n', res(2, i));
    fprintf('POD (L1): %0.4f\n', res(3, i));
    fprintf('POD (L2): %0.4f\n', res(4, i));
    fprintf('Training (L2): \t%0.4f\n', res(5, i));
    disp('---------------------------')
end
            
% Calculate error bars on the residual (assumes all samples are independent)
var_res = res2 - res.^2;  % Variance
err = sqrt(var_res/(num_configs*size(Test, 2)));

save('../output/re100-cylinder/reconst_loop_out.mat', 'res', 'err', 'sigma', 'num_configs', 'mTest')
