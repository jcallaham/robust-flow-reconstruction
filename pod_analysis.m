%% Calculate POD spectra and subspace projections for all data sets
% (mean subtracted)
%
% Jared Callaham 2018

clear all; close all; clc;
addpath('utils')

%% Mixing layer: linearly growing windows as a function of downstream coordinate
% Same windows are used here as in figure 9
disp('Mixing Layer: Windows')
disp('=====================================')

load('data/mixing_vort.mat', 'vort')
[ny, nx, m] = size(vort);
vort = reshape(vort, [], m); 

% Build linearly growing windows
window_width = @(x) floor(0.5*(x/nx)*ny);  % Linearly growing window width

x0 = round(0.14*nx);  % Beginning of first window
window_end = x0 + window_width(x0);  % Keep track of where the last window is
while window_end < nx  % Note that last one will be too far (stop the reconstruction loop short)
    x0 = [x0 window_end]; 
    window_end = window_end + window_width(window_end);
end
k = length(x0)-1;  % Number of kernels/windows
ns = 10;  % Evenly spaced along midline

% Partition data set
mTrain = round(0.9*size(vort, 2));
Train = vort(:, 1:mTrain);
Test = vort(:, mTrain+1:end);

m_mix = 2:10:size(Train,2);  % Variable number of training snapshots to compare
window_res = zeros(length(m_mix), k);   % Residual error as a function of m_mix
for i=1:length(m_mix)
    % Loop through windows
    for j=1:k
        window =  ny*x0(j):ny*(x0(j+1))-1;  % Indices of window in 2D matrix
        
        [U, ~, ~] = svd(Train(window, 1:m_mix(i)), 'econ');  % Compute POD basis
        res = Test(window, :) - U*(U'*Test(window, :));   % Residual errors
        res_denom = sum( Test(window, :).^2, 1);  % For normalizing residual errors
        
        window_res(i, j) = mean(sum(res.^2, 1)./res_denom);  % Normalized L2 error (3.10)
        fprintf('m: %d\t window: %d/%d\t res: %0.4f\n', m_mix(i), j, k, window_res(i, j));
    end
end
clearvars vort Test Train


%% Mixing Layer
disp('Mixing Layer')
disp('=====================================')

load('data/mixing_vort.mat', 'vort')
vort = reshape(vort, [], size(vort, 3));

% Calculate POD eigenvalues of full data set
[~, Sigma, ~] = svd(bsxfun(@minus, vort, mean(vort, 2)), 'econ');

mix_sv = diag(Sigma);  % Save singular values

% Partition data set
mTrain = round(0.9*size(vort, 2));
Train = vort(:, 1:mTrain);
Test = vort(:, mTrain+1:end);

n = size(Train, 1);

% Calculate subspace residual
m_mix = 2:10:size(Train,2);
mix_res = zeros(length(m_mix), 1);
rand_mix = zeros(length(m_mix), 1);
for i=1:length(m_mix)
    % Calculate average of POD subspace projection residual over test data
    mean_flow = mean(Train(:, 1:m_mix(i)), 2);  % Mean flow (given what we know so far)
    train_sub = bsxfun(@minus, Train(:, 1:m_mix(i)), mean_flow); % Mean-subtracted training set
    test_sub = bsxfun(@minus, Test, mean_flow); % Mean-subtracted test set
    
    [U, ~, ~] = svd(train_sub, 'econ');  % POD basis from training data
    res = test_sub - U*(U'*test_sub);  % Residual errors
    res_denom = sum( test_sub.^2, 1);  % For normalizing residual errorss
    mix_res(i) = mean(sum(res.^2, 1)./res_denom);

    % Random matrices of the same size for comparison
    [Q, ~] = qr(randn(n, m_mix(i)), 0);
    
    res = test_sub - Q*(Q'*test_sub);
    res = sum(res.^2, 1)./res_denom;
    rand_mix(i) = mean(res);
    
    fprintf('m: %d\t res: %0.4f\t rand: %0.4f\n', m_mix(i), mix_res(i), rand_mix(i));
end
clearvars vort Test Train

%% Sea Surface Temperature
disp('SST')
disp('=====================================')
load('data/sst_weekly.mat', 'sst')
mask = any(isnan(sst), 2);
sst(mask, :) = 0;
n = size(sst, 1);

% Calculate POD eigenvalues of full mean-subtracted data set
[~, Sigma, ~] = svd(bsxfun(@minus, sst, mean(sst, 2)), 'econ');
sst_sv = diag(Sigma);  % Singular values

% Partition data set
mTrain = 20*52;  % 20 years of training data
Train = sst(:, 1:mTrain);
Test = sst(:, mTrain+1:end);

m_sst = 2:10:size(Train,2);  % Number of training snapshots to use for subspace calculation
sst_res = zeros(length(m_sst), 1);
rand_sst = zeros(length(m_sst), 1);
for i=1:length(m_sst)
    % Calculate average of POD subspace projection residual over test data
    mean_flow = mean(Train(:, 1:m_sst(i)), 2);  % Mean flow (given what we know so far)
    train_sub = bsxfun(@minus, Train(:, 1:m_sst(i)), mean_flow); % Mean-subtracted training set
    test_sub = bsxfun(@minus, Test, mean_flow); % Mean-subtracted test set
    
    [U, ~, ~] = svd(train_sub, 'econ');  % Compute POD basis
    res = test_sub - U*(U'*test_sub);  % Residual errors
    res_denom = sum( test_sub.^2, 1);  % For normalizing residual errors... could also define using Test
    sst_res(i) = mean(sum(res.^2, 1)./res_denom);

    % Random matrices of the same size for comparison
    [Q, ~] = qr(randn(n, m_sst(i)), 0);
    
    res = test_sub - Q*(Q'*test_sub);
    res = sum(res.^2, 1)./res_denom;
    rand_sst(i) = mean(res);
        
    fprintf('m: %d\t res: %0.4f\t rand: %0.4f\n', m_sst(i), sst_res(i), rand_sst(i));
end

clearvars sst Test Train



%% Re=100 cylinder
disp('Cylinder')
disp('=====================================')
load('data/cylinder_vort.mat', 'VORTALL')

% Calculate POD eigenvalues of full mean-subtracted data set
[~, Sigma, ~] = svd(bsxfun(@minus, VORTALL, mean(VORTALL, 2)), 'econ');
cyl_sv = diag(Sigma);  % Singular values (to save to output)
n = size(VORTALL, 1);

% Partition data set (no mean subtraction)
mTrain = round(0.8*size(VORTALL, 2));
Train = VORTALL(:, 1:mTrain);
Test = VORTALL(:, mTrain+1:end); 

m_cyl = 2:mTrain;  % Number of snapshots to use in training
cyl_res = zeros(length(m_cyl), 1);
rand_cyl = zeros(length(m_cyl), 1);
for i=1:length(m_cyl)    
    % Calculate average of POD subspace projection residual over test data
    mean_flow = mean(Train(:, 1:m_cyl(i)), 2);  % Mean flow (given what we know so far)
    train_sub = bsxfun(@minus, Train(:, 1:m_cyl(i)), mean_flow); % Mean-subtracted training set
    test_sub = bsxfun(@minus, Test, mean_flow); % Mean-subtracted test set
    
    [U, ~, ~] = svd(train_sub, 'econ');  % Compute POD basis
    res = test_sub - U*(U'*test_sub);  % Residual errors
    res_denom = sum( test_sub.^2, 1);  % For normalizing residual errors... could also define using Test
    cyl_res(i) = mean(sum(res.^2, 1)./res_denom);

    % Random matrices of the same size for comparison
    [Q, ~] = qr(randn(n, m_cyl(i)), 0);
    
    res = test_sub - Q*(Q'*test_sub);
    res = sum(res.^2, 1)./res_denom;
    rand_cyl(i) = mean(res);
    
    fprintf('m: %d\t res: %0.4f\t rand: %0.4f\n', m_cyl(i), cyl_res(i), rand_cyl(i));
end

clearvars VORTALL Test Train


%% HYCOM Gulf of Mexico vorticity (and comparison to random)
disp('HYCOM Gulf of Mexico')
disp('=====================================')

load('data/gom_vort.mat', 'vort')

% Calculate POD eigenvalues of full data set (without NaN values)
mask = any(isnan(vort), 2);
vort(mask, :) = 0;
[~, Sigma, ~] = svd(bsxfun(@minus, vort, mean(vort, 2)), 'econ');

gom_sv = diag(Sigma);  % Singular values (to save)

% Partition data set and calculate subspace projection residual
mTrain = round(0.9*size(vort, 2));
Train = vort(:, 1:mTrain);
Test = vort(:, mTrain+1:end);

n = size(Train, 1);

m_gom = 2:10:size(Train,2);
gom_res = zeros(length(m_gom), 1);
rand_gom = zeros(length(m_gom), 1);
res_denom = sum( Test.^2, 1);  % For normalizing residual errors
for i=1:length(m_gom)
    % Calculate average of POD subspace projection residual over test data
    mean_flow = mean(Train(:, 1:m_gom(i)), 2);  % Mean flow (given what we know so far)
    train_sub = bsxfun(@minus, Train(:, 1:m_gom(i)), mean_flow); % Mean-subtracted training set
    test_sub = bsxfun(@minus, Test, mean_flow); % Mean-subtracted test set
    
    [U, ~, ~] = svd(train_sub, 'econ');  % Compute POD basis
    res = test_sub - U*(U'*test_sub);  % Residual errors
    res_denom = sum( test_sub.^2, 1);  % For normalizing residual errors... could also define using Test
    gom_res(i) = mean(sum(res.^2, 1)./res_denom);

    % Random matrices of the same size for comparison
    [Q, ~] = qr(randn(n, m_gom(i)), 0);
    
    res = test_sub - Q*(Q'*test_sub);
    res = sum(res.^2, 1)./res_denom;
    rand_gom(i) = mean(res);
    
    fprintf('m: %d\t res: %0.4f\t rand: %0.4f\n', m_gom(i), gom_res(i), rand_gom(i));
end
clearvars vort Test Train


%% Save results
save('output/sv_spectra.mat', 'cyl_sv', 'sst_sv', 'gom_sv', 'mix_sv')

save('output/subspace_residuals.mat', 'cyl_res', 'sst_res', 'gom_res'...
    , 'm_cyl', 'm_sst', 'm_gom', 'rand_cyl', 'rand_sst', 'rand_gom'...
    , 'm_mix', 'mix_res', 'rand_mix', 'window_res')
disp('Done')

