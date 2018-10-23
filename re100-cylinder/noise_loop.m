%% Error vs Gaussian noise level for different measurement strategies
%  (all using L1 optimization and training dictionary)
%
% Jared Callaham, 2018
warning('off', 'all');
addpath('../utils')

load_cylinder;  % Load full data set and dimensions
define_measurements;  % Define or load measurement matrices C

rms_vort = flow.avg_energy;
mTrain = size(Train, 2);

energy_rescale = true;

%% Define dictionaries (just for efficiency)
D_window = C_window*Train; 
D_rand = C_rand*Train;
D_vslice = C_vslice*Train;
D_hslice = C_hslice*Train;

%% Loop through corruption levels
sigma = 0:0.01:0.2;     % Gaussian noise level

res = zeros(4, length(sigma));  % Normalized residual error
res2 = zeros(4, length(sigma));  % Same, but squared (to calculate variance)
num_iters = 10;   % Number of noise realizations to test (actually have mTest*num_iters)

disp('Beginning loop...')
for i=1:length(sigma)
    for j=1:num_iters
        fprintf('%d ', j);
        for t=1:size(Test, 2)
            x = Test(:, t);
            % Corrupt with noise
            noise = sigma(i)*rms_vort*randn(size(x));
            
            %% Optimized reconstruction (training dictionary, L1 norm)
            % Window
            y = C_window*(x + noise); % Measure noisy slice
            % Sparse representation in terms of the dictionary
            s = sp_approx(y, D_window, sigma(i), flow);  
            % Residual error in reconstruction (use only coefficients, not identified noise)            
            [~, r] = reconstruct(x, Train, s, flow, energy_rescale);
            res(1, i) = res(1, i) + r;
            res2(1, i) = res2(1, i) + r^2;

            % Random points
            y = C_rand*(x + noise); % Noisy slice
            
            % Sparse representation in terms of the dictionary
            s = sp_approx(y, D_rand, sigma(i), flow);  
            % Residual error in reconstruction
            [~, r] = reconstruct(x, Train, s, flow, energy_rescale);
            res(2, i) = res(2, i) + r;
            res2(2, i) = res2(2, i) + r^2;
            
            % Horizontal slice
            y = C_hslice*(x + noise); % Noisy slice
            
            % Sparse representation in terms of the dictionary
            s = sp_approx(y, D_hslice, sigma(i), flow);  
            % Residual error in reconstruction
            [~, r] = reconstruct(x, Train, s, flow, energy_rescale);
            res(3, i) = res(3, i) + r;
            res2(3, i) = res2(3, i) + r^2;
            
            % Vertical slice
            y = C_vslice*(x + noise); % Noisy slice
            % Sparse representation in terms of the dictionary
            s = sp_approx(y, D_vslice, sigma(i), flow);  
            % Residual error in reconstruction
            [~, r] = reconstruct(x, Train, s, flow, energy_rescale);
            res(4, i) = res(4, i) + r;
            res2(4, i) = res2(4, i) + r^2;
        end
    end
    fprintf('\n');
    res(1:4, i) = res(1:4, i)/(num_iters*size(Test, 2));    % Avg error
    res(5, i) = res(5, i)/(num_iters*(size(Test, 2)-delays+1)); % Adjust separately for delayed case    
    res2(1:4, i) = res2(1:4, i)/(num_iters*size(Test, 2));  % Avg error, squared 
    fprintf('Noise level: \t%0.2f\n', sigma(i))
    fprintf('Window: \t%0.4f\n', res(1, i));
    fprintf('Rand pts: \t%0.4f\n', res(2, i));
    fprintf('Horz slice: %0.4f\n', res(3, i));
    fprintf('Vert slice: \t%0.4f\n', res(4, i));
end
            
% Calculate error bars on the residual (assumes all samples are independent)
var_res = res2 - res.^2;  % Variance
err = sqrt(var_res/(num_iters*size(Test, 2)));

save('../output/re100-cylinder/noise_loop_out.mat', 'res', 'err', 'sigma', 'num_iters')
