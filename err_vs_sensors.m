%% Compare reconstruction error vs number of random point measurements for all data sets
% (global reconstructions, noise added for cylinder and mixing layer)
%
% Jared Callaham 2018

clear all; close all; clc;
addpath(genpath('utils'))

ns = [1 3 6 round(logspace(1, 4, 20))];  % Number of point measurements
sigma = [0, 0.1, 0.2, 0.3];   % Noise levels for cylinder and mixing layer (not in paper)
num_iters = 1;  % Number of independent realizations of sensor placement



%% NOAA SST
disp('Sea surface temperature')
disp('=========================================================');
run('noaa-sst/load_sst.m');
ny = flow.ny; mask = flow.mask;
Train(mask, :) = 0; Test(mask, :) = 0;

% Loop through number of measurements
K = 5; % OMP sparsity (estimated empirically)
sst_res = zeros([length(ns), 1]);
res2 = zeros([length(ns), 1]);
for ns_idx=1:length(ns)
    for test_idx=1:size(Test, 2)
        for iter=1:num_iters
            sensor_idx = sst_sensors(ns(ns_idx), mask, ny, lat);

            D = double(Train(sensor_idx, :));  % Construct dictionary

            % Normalize library columns for OMP
            omp_rescale = zeros(size(Train, 2), 1);
            for j=1:size(Train, 2)
                omp_rescale(j) = norm(D(:, j));
                D(:, j) = D(:, j)/omp_rescale(j);
            end

            x = double(Test(:, test_idx));      % Test field
            y = x(sensor_idx);  % Measure field
            s = omp(D, y, D'*D, K); % Sparse approximation via OMP

            [~, r] = reconstruct(x, Train, full(s)./omp_rescale, flow, false); 
            if ~isnan(r)
                sst_res(ns_idx) = sst_res(ns_idx) + r;
                res2(ns_idx) = res2(ns_idx) + r^2;
            end

        end
    end
    % Average and display status
    sst_res(ns_idx) = sst_res(ns_idx)/(size(Test, 2)*num_iters);
    res2(ns_idx) = res2(ns_idx)/(size(Test, 2)*num_iters);
    fprintf('ns: %i\t res: %0.4f\n', ns(ns_idx), sst_res(ns_idx));
end

% Calculate error
var_res = res2' - sst_res.^2;  % Variance
sst_err = sqrt(var_res/(num_iters*size(Test, 2)));

%% HYCOM GoM (no noise)
disp('HYCOM Gulf of Mexico')
disp('=========================================================');
% Load and partition data
run('hycom/load_hycom.m');
nx = flow.nx; ny = flow.ny; mask = flow.mask;
clearvars vort

% Loop through number of measurements
K = 400; % OMP sparsity (very rough empirical estimate)
gom_res = zeros([length(ns), 1]);
res2 = zeros([length(ns), 1]);
for ns_idx=1:length(ns)
    for test_idx=1:size(Test, 2)
        for iter=1:num_iters 
            % Choose randomly from valid locations
            sensor_idx = randsample(find(~mask), ns(ns_idx)); 
            
            x = double(Test(:, test_idx));      % Test field
            y = x(sensor_idx);  % Measure field
            
            D = double(Train(sensor_idx, :));  % Construct dictionary

            % Normalize library columns for OMP
            omp_rescale = zeros(size(Train, 2), 1);
            for j=1:size(Train, 2)
                omp_rescale(j) = norm(D(:, j));
                D(:, j) = D(:, j)/omp_rescale(j);
            end
            
            s = omp(D, y, D'*D, K); % Sparse approximation via OMP

            [~, r] = reconstruct(x, Train, full(s)./omp_rescale, flow, false); 
            if ~isnan(r)
                gom_res(ns_idx) = gom_res(ns_idx) + r;
                res2(ns_idx) = res2(ns_idx) + r^2;
            end

        end
    end
    % Average and display status
    gom_res(ns_idx) = gom_res(ns_idx)/(size(Test, 2)*num_iters);
    res2(ns_idx) = res2(ns_idx)/(size(Test, 2)*num_iters);
    fprintf('ns: %i\t res: %0.4f\n', ns(ns_idx), gom_res(ns_idx));
end

% Calculate standard error
var_res = res2' - gom_res.^2;  % Variance
gom_err = sqrt(var_res/(num_iters*size(Test, 2)));



%% Re=100 cylinder
disp('Re=100 cylinder')
disp('=========================================================');
% Load and partition data
run('re100-cylinder/load_cylinder.m');
nx = flow.nx; ny = flow.ny; rms_vort = flow.avg_energy;
clearvars VORTALL

K = 1;  % OMP sparsity (here this is based on knowledge of periodicity)
cyl_res = zeros(length(sigma), length(ns));
res2 = zeros([length(ns), 1]);
cyl_err = zeros(size(cyl_res));
for noise_idx=1:length(sigma)
    for ns_idx=1:length(ns)
        for test_idx=1:size(Test, 2)
            for iter=1:num_iters
                sensor_idx = cylinder_sensors(ns(ns_idx), nx, ny);

                D = Train(sensor_idx, :);  % Construct dictionary

                % Normalize library columns for OMP
                omp_rescale = zeros(size(Train, 2), 1);
                for j=1:size(Train, 2)
                    omp_rescale(j) = norm(D(:, j));
                    D(:, j) = D(:, j)/omp_rescale(j);
                end

                x = Test(:, test_idx);      % Test field
                noise = sigma(noise_idx)*rms_vort*randn([ns(ns_idx), 1]);
                y = x(sensor_idx) + noise;  % Measure field
                s = omp(D, y, D'*D, K); % Sparse approximation via OMP

                [~, r] = reconstruct(x, Train, s./omp_rescale, flow, true);
                if ~isnan(r)
                    cyl_res(noise_idx, ns_idx) = cyl_res(noise_idx, ns_idx) + r;
                    res2(ns_idx) = res2(ns_idx) + r^2;
                end

            end
        end
        % Average and display status
        cyl_res(noise_idx, ns_idx) = cyl_res(noise_idx, ns_idx)/(size(Test, 2)*num_iters);
        res2(ns_idx) = res2(ns_idx)/(size(Test, 2)*num_iters);
        fprintf('sigma:%0.2f  ns: %i\t res: %0.4f\n'...
            , sigma(noise_idx), ns(ns_idx), cyl_res(noise_idx, ns_idx));
    end

    % Calculate error
    var_res = res2' - cyl_res(noise_idx, :).^2;  % Variance
    cyl_err(noise_idx, :) = sqrt(var_res/(num_iters*size(Test, 2)));
end

%% Mixing layer
disp('Mixing layer')
disp('=========================================================');
% Load and partition
run('mixing-layer/load_ml.m');
clearvars vort
nx = flow.nx; ny = flow.ny;
rms_vort = double(mean( sqrt( mean(Train.^2, 1) ) ));  % For choosing noise levels

% Loop through number of measurements
K = 72; % OMP sparsity (estimated empirically)
ml_res = zeros(length(sigma), length(ns));
res2 = zeros([length(ns), 1]);
ml_err = zeros(size(ml_res));
for noise_idx=1:length(sigma)
    for ns_idx=1:length(ns)
        for test_idx=1:size(Test, 2)
            for iter=1:num_iters
                % Choose randomly from valid locations
                sensor_idx = ml_sensors(ns(ns_idx), x_cc, y_cc);  

                D = double(Train(sensor_idx, :));  % Construct dictionary

                % Normalize library columns for OMP
                omp_rescale = zeros(size(Train, 2), 1);
                for j=1:size(Train, 2)
                    omp_rescale(j) = norm(D(:, j));
                    D(:, j) = D(:, j)/omp_rescale(j);
                end

                x = double(Test(:, test_idx));      % Test field
                noise = sigma(noise_idx)*rms_vort*randn([ns(ns_idx) 1]);
                y = x(sensor_idx) + noise;  % Measure field
                s = omp(D, y, D'*D, K); % Sparse approximation via OMP'
                [~, r] = reconstruct(x, Train, full(s)./omp_rescale, flow, true);
                if ~isnan(r)
                    ml_res(noise_idx, ns_idx) = ml_res(noise_idx, ns_idx) + r;
                    res2(ns_idx) = res2(ns_idx) + r^2;
                end

            end
        end
        % Average and display status
        ml_res(noise_idx, ns_idx) = ml_res(noise_idx, ns_idx)/(size(Test, 2)*num_iters);
        res2(ns_idx) = res2(ns_idx)/(size(Test, 2)*num_iters);
        fprintf('sigma:%0.2f  ns: %i\t res: %0.4f\n'...
            , sigma(noise_idx), ns(ns_idx), ml_res(noise_idx, ns_idx));
    end
    
    % Calculate error
    var_res = res2' - ml_res(noise_idx, :).^2;  % Variance
    ml_err(noise_idx, :) = sqrt(var_res/(num_iters*size(Test, 2)));
end


%% Save results
save('output/err_vs_ns.mat', 'ns', 'noise', 'cyl_res', 'cyl_err'...
    , 'ml_res', 'ml_err', 'sst_res', 'sst_err', 'gom_res', 'gom_err');



%% ========================================================
% Functions to choose random (but reasonable) sensor locations

function sensor_idx = ml_sensors(ns, x_cc, y_cc)
    % Chooses random points in linear growth of mixing layer
    sensor_idx = zeros(ns, 1);
    for i=1:ns
        found = false;
        while(~found)
            temp_idx = randi(length(x_cc)*length(y_cc)-1);
            yloc = 1+mod(temp_idx, length(y_cc)); % y as grid location
            xloc = 1+floor(temp_idx/length(y_cc));    % x as grid location
            yloc = y_cc(yloc); xloc = x_cc(xloc);   % convert to flow coordinates
            found = ~ismember(temp_idx, sensor_idx) && (160*abs(yloc) < 6*xloc);  % Check that location is inside layer
        end
        sensor_idx(i) = temp_idx;
    end
end

function sensor_idx = cylinder_sensors(ns, nx, ny)
    % Chooses random points in wake of cylinder (back 80% and middle 50%)
    sensor_idx = sort(randperm(floor(0.8*ny)*floor(0.5*nx), ns));
    sensor_idx = [floor(sensor_idx/floor(0.5*nx))' mod(sensor_idx, floor(0.5*nx))']; % Convert to 2D locations
    sensor_idx = [floor(0.2*ny)+sensor_idx(:, 1) floor(0.25*nx)+sensor_idx(:, 2)];  % Translate to wake
    sensor_idx = sensor_idx(:, 2) + (nx-1)*sensor_idx(:, 1);  % Convert back to indices
end

function sensor_idx = sst_sensors( ns, mask, ny, lat )
    % Choose random valid measurement locations on (50S, 50N)
    n = length(mask);
    sensor_idx = zeros(ns, 1);
    for i=1:ns
        found = false;
        while(~found)
            temp_idx = randi(n);
            yloc = lat(mod(temp_idx, ny)+1);  % Latitude corresponding to y-coordinat in grid
            found = ~mask(temp_idx) && (yloc > -50) && (yloc < 50)...
                && ~ismember(temp_idx, sensor_idx);
        end
        sensor_idx(i) = temp_idx;
    end

end
