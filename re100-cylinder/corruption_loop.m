%% Error vs sparse corruption for different measurement strategies
%  (all using L1 optimization and training dictionary)
%
% Jared Callaham, 2018
warning('off', 'all');
addpath('../utils')

load_cylinder;  % Load and partition full data set
define_measurements;  % Define or load measurement matrices C

mTrain = size(Train, 2);
rms_vort = flow.avg_energy;

energy_rescale = true;


%% Define dictionaries (just for efficiency)
% Extended dictionary to handle corruption
ns = size(C_window, 1);
D_window = [C_window*Train, spdiags(ones(ns, 1), 0, ns, ns)]; 
ns = size(C_rand, 1);
D_rand = [C_rand*Train, spdiags(ones(ns, 1), 0, ns, ns)]; 
ns = size(C_vslice, 1);
D_vslice = [C_vslice*Train, spdiags(ones(ns, 1), 0, ns, ns)]; 
ns = size(C_hslice, 1);
D_hslice = [C_hslice*Train, spdiags(ones(ns, 1), 0, ns, ns)]; 

%% Loop through corruption levels
sigma = 0.01;     % Gaussian noise level
rho = 0:0.05:1;   % Percentage of corrupt measurements
max_noise = max(max(abs(bsxfun(@plus, Train, mean_flow))));

res = zeros(4, length(rho));  % Normalized residual error
res2 = zeros(4, length(rho));  % Same, but squared (to calculate variance)
num_iters = 10;   % Number of noise realizations to test (actually have mTest*num_iters)

disp('Beginning loop...')
for i=1:length(rho)
    for j=1:num_iters
        fprintf('%d ', j);
        for t=1:size(Test, 2)
            x = Test(:, t);
            % Corrupt with noise
            noise = sigma*rms_vort*randn(size(x));
            corrupt_idx = randperm(n, round(rho(i)*n));  % Choose random locations to corrupt
            corruption = 2*max_noise*(0.5 - rand(size(corrupt_idx)));
            x_corrupt = x;
            x_corrupt(corrupt_idx) = corruption;  % Replace corrupted locations with noise
            
            %% Optimized reconstruction (training dictionary, L1 norm)
            % Window
            y = C_window*(x_corrupt + noise); % Measure noisy slice
            % Sparse representation in terms of the dictionary
            w = sp_approx(y, D_window, sigma, flow);  
            % Residual error in reconstruction (use only coefficients, not identified noise)            
            [~, r] = reconstruct(x, Train, w(1:mTrain), flow, energy_rescale);
            res(1, i) = res(1, i) + r;
            res2(1, i) = res2(1, i) + r^2;

            % Random points
            y = C_rand*(x_corrupt + noise); % Noisy slice
            
            % Sparse representation in terms of the dictionary
            w = sp_approx(y, D_rand, sigma, flow);  
            % Residual error in reconstruction
            [~, r] = reconstruct(x, Train, w(1:mTrain), flow, energy_rescale);
            res(2, i) = res(2, i) + r;
            res2(2, i) = res2(2, i) + r^2;
            
            % Horizontal slice
            y = C_hslice*(x_corrupt + noise); % Noisy slice
            
            % Sparse representation in terms of the dictionary
            w = sp_approx(y, D_hslice, sigma, flow);  
            % Residual error in reconstruction
            [~, r] = reconstruct(x, Train, w(1:mTrain), flow, energy_rescale);
            res(3, i) = res(3, i) + r;
            res2(3, i) = res2(3, i) + r^2;
            
            % Vertical slice
            y = C_vslice*(x_corrupt + noise); % Noisy slice
            % Sparse representation in terms of the dictionary
            w = sp_approx(y, D_vslice, sigma, flow);  
            % Residual error in reconstruction
            [~, r] = reconstruct(x, Train, w(1:mTrain), flow, energy_rescale);
            res(4, i) = res(4, i) + r;
            res2(4, i) = res2(4, i) + r^2;
        end
        
    end
    fprintf('\n');
    res(1:4, i) = res(1:4, i)/(num_iters*size(Test, 2));    % Avg error
    res2(1:4, i) = res2(1:4, i)/(num_iters*size(Test, 2));  % Avg error, squared 
    fprintf('Corruption level: \t%0.2f\n', rho(i))
    fprintf('Window: \t%0.4f\n', res(1, i));
    fprintf('Rand pts: \t%0.4f\n', res(2, i));
    fprintf('Horz slice: %0.4f\n', res(3, i));
    fprintf('Vert slice: \t%0.4f\n', res(4, i));
end
            
% Calculate error bars on the residual (assumes all samples are independent)
var_res = res2 - res.^2;  % Variance
err = sqrt(var_res/(num_iters*size(Test, 2)));

save('../output/re100-cylinder/corrupt_loop_out.mat', 'res', 'err', 'sigma', 'rho', 'num_iters')
