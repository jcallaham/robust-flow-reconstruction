%% Load cylinder data files and useful parameters
% Jared Callaham 2018

load('../data/cylinder_vort.mat')

flow.nx = nx; flow.ny = ny;
flow.n = nx*ny;  % y is streamwise direction, x is cross-stream
flow.m = size(VORTALL, 2);  % Total number of available snapshots

flow.mTrain = 32;  % Number of snapshots to use for training (one full period of vortex shedding)

% Partition data into training and test set. Also returns the mean flow and any masked values
[Train, Test, flow.mean_flow, ~] = partition(VORTALL, flow.mTrain);

% RMS vorticity in the training set - used for rescaling
flow.avg_energy = mean( sqrt( mean(Train.^2, 1) ) );

% Helpful for plotting
load('../utils/CCcool.mat')
flow.cmap = CC;
flow.clim = [-5, 5];