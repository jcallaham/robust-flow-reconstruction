%% Load HYCOM Gulf of Mexico velocity data and flow parameters
% Jared Callaham 2018


load('../data/gom_vort.mat')

flow.nx = length(lon); flow.ny = length(lat);
flow.n = flow.nx*flow.ny;  %y is streamwise direction, x is spanwise
flow.m = size(vort, 2);

[Train, Test, flow.mean_flow, flow.mask] = partition(vort, 0.9);
flow.mean_flow = double(flow.mean_flow);
flow.avg_energy = NaN;  % No energy rescaling
flow.mTrain = size(Train, 2);

% Helpful for plotting
load('../utils/CCcool.mat')
flow.cmap = CC;
flow.clim = [-10, 10];

clearvars vort