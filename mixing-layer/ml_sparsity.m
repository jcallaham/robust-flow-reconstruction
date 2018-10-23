%% Compare sparsity of mixing layer reconstruction with downstream windowing to global reconstruction
% Centerline measurements

disp('Initializing... ')
addpath(genpath('../utils'))
load_ml;  % Load and partition mixing layer
nx = flow.nx; ny = flow.ny;
[XX, YY] = meshgrid(x_cc, y_cc);

%% Construct linearly growing windows as a function of downstream coordinate

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

% Find window and global sensor locations
sensor_idx = []; % All sensor locations (for later global reconstruction)
for i=1:k
    % Stagger sensors along midline of window
    x_span = floor( (1:ns)*(x0(i+1)-x0(i))/(ns+1) );  % Evenly space sensors along midline of window
    sensor_locs = [x_span + x0(i); ny/2+0*x_span]';   % Horizontal slice
    sensor_idx = [sensor_idx (sensor_locs(:, 1)*ny + sensor_locs(:, 2))'];  % Save sensor indices for global reconstruction
end


%% Reconstruction
% Sparsity of representation
t = 3:2:size(Test, 2);
K_global_sr = zeros([length(t), 1]);
K_window_sr = zeros([length(t), 1]);

% Residual reconstruction error
res_global_sr = zeros([length(t), 1]);
res_window_sr = zeros([length(t), 1]);

disp('Beginning loop...')
for i=1:length(t) % Choose test snapshot
    x = double(Test(:, t(i)));
    
    sigma = 0.3;  % Noise level
    noise = sigma*flow.avg_energy*randn([flow.n, 1]);
    
    % Global estimations
    x_window_sr = zeros(size(x));
    for j =1:k   % Loop through windows

        %% Instantaneous window reconstruction
        window =  ny*x0(j):ny*(x0(j+1))-1;  % Window for comparison
        window_sensors = sensor_idx(ismember(sensor_idx, window)); % Find sensors in window
        
        local_flow = flow;
        local_flow.avg_energy = mean( sqrt( mean(Train(window, :).^2, 1) ) );
        
        D = double(Train(window_sensors, :));
        y = x(window_sensors) + noise(window_sensors);
        
        s = sp_approx(y, D, sigma, flow);
        [~, r] = reconstruct(x(window), Train(window, :), s, local_flow, true);

        % Update relative sparsity
        K_window_sr(i) = K_window_sr(i) + sum(s>1e-6)/(flow.mTrain*k); 

    end
    res_window_sr(i) = norm(x_window_sr - x)/norm(x);

    %% Global reconstruction
    D = double(Train(sensor_idx, :));
    y = x(sensor_idx) + noise(sensor_idx);

    s = sp_approx(y, D, sigma, flow);
    % Reconstruct (zeros for mean, since this wasn't subtracted)
    [~, res_global_sr(i)] = reconstruct(x, Train, s, flow, true);
    K_global_sr(i) = sum(abs(s)>1e-6)/flow.mTrain;
    
    
    fprintf('%d/%d\n', i, length(t))
end

%% Save results

save('../output/mixing-layer/sparsity.mat', 'flow', 'sigma', 'K_window_sr'...
    , 'K_global_sr', 'res_window_sr', 'res_global_sr', 'k')
    
