%% Load cylinder data files and useful parameters
% Jared Callaham 2018

load('../data/sst_weekly.mat')

flow.nx = length(lon); flow.ny = length(lat);
flow.n = flow.nx*flow.ny;  %y is streamwise direction, x is spanwise
flow.m = size(sst, 2);

flow.mTrain = 20*52;  % 20 years of weekly training data

[Train, Test, ~, flow.mask] = partition(sst, flow.mTrain);
flow.mean_flow = 0;  % Mean is subtracted, but we won't add it back in

% Calculate rms temperature fluctuations
flow.avg_energy = double(mean( sqrt( mean(Train.^2, 1) ) ));

% Helpful for plotting
load('../utils/CCcool.mat')
flow.cmap = CC;
flow.clim = [-10, 10];

clearvars sst