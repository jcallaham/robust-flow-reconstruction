%% Compare sparsity of mixing layer reconstruction with downstream windowing to global reconstruction


disp('Initializing... ')
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

% Find sensor locations
sensor_idx = []; % All sensor locations (for later global reconstruction)
for i=1:k
    % Stagger sensors along midline of window
    x_span = floor( (0:ns-1)*window_width(x0(i))/ns );  % Evenly space sensors along midline of window
    sensor_locs = [x_span + x0(i); ny/2+0*x_span]';   % Horizontal slice
    sensor_idx = [sensor_idx (sensor_locs(:, 1)*ny + sensor_locs(:, 2))'];  % Save sensor indices for global reconstruction
end

%% POD library
disp('Computing POD... ')
[U, Sigma, ~] = svd(Train, 'econ');

beta = 0.5*size(Train,2)/size(Train,1); % Aspect ratio of data matrix
thresh = optimal_SVHT_coef(beta,0) * median(diag(Sigma));
r_GD = find(diag(Sigma)>thresh, 1, 'last');
Ur = double(U(:, 1:r_GD));  % Truncate
clearvars U % Free up memory
fprintf('Done\n')

%% Reconstruction
% Sparsity of representation
K_global_sr = zeros([size(Test, 2), 1]);
K_window_sr = zeros([size(Test, 2), 1]);
K_global_pod = zeros([size(Test, 2), 1]);
K_window_pod = zeros([size(Test, 2), 1]);
% Residual reconstruction error
res_global_sr = zeros([size(Test, 2), 1]);
res_window_sr = zeros([size(Test, 2), 1]);
res_global_pod = zeros([size(Test, 2), 1]);
res_window_pod = zeros([size(Test, 2), 1]);

% Total area covered by windows (for normalizing residual errors)
all_windowed = ny*x0(1):ny*(x0(end))-1;
windowed_points = length(all_windowed); 

disp('Beginning loop...')
for t=1:size(Test, 2) % Choose test snapshot
    x = double(Test(:, t));
    
    sigma = 0.7;  % Noise level
    rms_vort = double(mean( sqrt( mean(Train.^2, 1) ) ));  % For noise level only - no energy rescaling
    noise = sigma*rms_vort*randn([flow.n, 1]);
    
    % Global estimations
    x_window_sr = zeros(size(x));
    x_window_pod = zeros(size(x));
    for j =1:k   % Loop through windows

        %% Instantaneous window reconstruction
        window =  ny*x0(j):ny*(x0(j+1))-1;  % Window for comparison
        window_sensors = sensor_idx(ismember(sensor_idx, window)); % Find sensors in window
        D = double(Train(window_sensors, :));
        y = x(window_sensors) + noise(window_sensors);
        eps = 3*sigma*rms_vort*sqrt(length(y));  % Tolerance for optimization problem

        s = sp_approx(y, D, sigma, flow);
        [~, r] = reconstruct(x(window), Train(window, :), s, flow);

        % Update relative sparsity
        K_window_sr(t) = K_window_sr(t) + sum(s>1e-6)/(flow.mTrain*k);
        % Update error (weighted by window size)
        res_window_sr(t) = res_window_sr(t) + r*length(window)/windowed_points;  

        %% L2 POD reconstruction in windows
        s = lsqminnorm(Ur(window_sensors, :), y);
        [~, r] = reconstruct(x(window), Ur(window, :), s, flow);
        
        % Update relative sparsity
        K_window_pod(t) = K_window_pod(t) + sum(abs(s)>1e-6)/(r_GD*k);
        % Update error (weighted by window size)
        res_window_pod(t) = res_window_pod(t) + r*length(window)/windowed_points;  
    end

    %% Global reconstruction
    D = double(Train(sensor_idx, :));
    y = x(sensor_idx) + noise(sensor_idx);

    s = sp_approx(y, D, sigma, flow);
    % Reconstruct (zeros for mean, since this wasn't subtracted)
    [x_hat, res_global_sr(t)] = reconstruct(x, Train, s, flow);
    K_global_sr(t) = sum(abs(s)>1e-6)/flow.mTrain;

    % L2 POD
    s = lsqminnorm(Ur(sensor_idx, :), y);
    [x_hat, res_global_pod(t)] = reconstruct(x, Ur, s, flow);
    K_global_pod(t) = sum(abs(s)>1e-6)/r_GD;
    
    
    fprintf('%d/%d\n', t, size(Test, 2))
    fprintf('Residuals\n\tWindow SR:%0.2f  Window POD:%0.2f  Global SR:%0.2f  Global POD:%0.2f\n'...
        , res_window_sr(t), res_window_pod(t), res_global_sr(t), res_global_pod(t));
    fprintf('Sparsity\n\tWindow SR:%0.4f  Window POD:%0.2f  Global SR:%0.4f  Global POD:%0.2f\n'...
        , K_window_sr(t), K_window_pod(t), K_global_sr(t), K_global_pod(t));

end

%% Save results

save('../output/mixing-layer/sparsity.mat', 'flow', 'sigma', 'K_window_sr'...
    , 'K_window_pod', 'K_global_sr', 'K_global_pod', 'res_window_sr'...
    , 'res_window_pod', 'res_global_sr', 'res_global_pod')
    
