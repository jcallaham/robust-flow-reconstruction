%% Load cylinder data files and useful parameters
% Jared Callaham 2018

load('../data/mixing_vort.mat')

[flow.ny, flow.nx, flow.m] = size(vort);

%% Downstream cutoff (optional)
flow.nx = floor(0.8*flow.nx);
vort = vort(:, 1:flow.nx, :);
x_cc = x_cc(1:flow.nx);

flow.n = flow.nx*flow.ny;  %y is streamwise direction, x is spanwise
vort = reshape(vort, flow.n, flow.m);  % Stack snapshots into columns

%% Partition (no mean subtraction or masking)
flow.mTrain = floor(0.9*flow.m);
Train = vort(:, 1:flow.mTrain);
Test = vort(:, flow.mTrain+1:end);
flow.mean_flow = 0;  % No mean subtraction
flow.avg_energy = double(mean( sqrt( mean(Train.^2, 1) ) ));

%% Helpful for plotting
load('../utils/cmap_marine.mat', 'marine_full')
flow.cmap = marine_full;
flow.clim = [-.5, .5];

clearvars vort