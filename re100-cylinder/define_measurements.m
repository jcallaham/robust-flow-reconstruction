%% Define measurement locations and matrices for different measurement strategies
% Jared Callaham, 2018
% load('../data/cylinder_vort.mat', 'nx', 'ny');
% n = nx*ny;
% 
% %% Define measurement
% % Random point measurements
% ns = 10;  % Number of point measurements
% 
% % Restrict to cylinder wake: final 80% of width, middle 50% of height
% %rand_loc = [randperm(round(0.8*ny), ns)' randperm(round(0.5*nx), ns)']; % Choose sensors on restricted area
% %rand_loc = [round(0.2*ny)+rand_loc(:, 1) round(0.25*nx)+rand_loc(:, 2)];  % Translate to wake
% C_rand = spdiags(ones(n, 1), 0, n, n);
% C_rand = C_rand(rand_loc(:, 2) + (nx-1)*rand_loc(:, 1), :);  % Sparse measurement matrix
% 
% % Horizontal slice
% hslice_loc = [(165:360)', 0*(165:360)' + 100];
% C_hslice = spdiags(ones(n, 1), 0, n, n);
% C_hslice = C_hslice(hslice_loc(:, 2) + (nx-1)*hslice_loc(:, 1), :);  % Measurement matrix
% 
% % Vertical slice
% vslice_loc = [0*(50:150)'+200, (50:150)'];
% C_vslice = spdiags(ones(n, 1), 0, n, n);
% C_vslice = C_vslice(vslice_loc(:, 2) + (nx-1)*vslice_loc(:, 1), :);  % Measurement matrix
% 
% % Window
% delta_x = round(0.3*nx);  % Spanwise window width
% delta_y = round(0.3*ny);  % Streamwise width
% ns = delta_x*delta_y;
% window_loc = 1:ns;
% window_loc = [round(window_loc/delta_x)' mod(window_loc, delta_x)' ]; % Convert to 2D locations
% window_loc = [round(0.4*ny)+window_loc(:, 1) round((nx-delta_x)/2)+window_loc(:, 2)];  % Translate to wake
% 
% % Convert to measurement matrix
% C_window = spdiags(ones(n, 1), 0, n, n);
% C_window = C_window(window_loc(:, 2) + (nx-1)*window_loc(:, 1), :);
% 
% clearvars ns nx ny n
% save('../output/re100-cylinder/measurement_matrices.mat')

load('../output/re100-cylinder/measurement_matrices.mat')
