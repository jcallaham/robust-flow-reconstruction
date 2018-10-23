% Run all scripts and generate figures
%
% Jared Callaham 2018

clear all; close all; clc;
addpath('utils')

%% Run scripts to generate results from raw data
% (note: it would take a really long time to run these all at once)
RUN_CODE = false;
if RUN_CODE
    % Compute singular value spectra (Fig 4) and subspace residuals (Fig 12)
    run('pod_analysis.m')  
    
    %% Re=100 cylinder
    % Measurement strategies and error vs noise (Fig 5)
    run('re100-cylinder/noise_loop.m')
    run('re100-cylinder/corruption_loop.m')
    
    % err vs noise for different dictionaries and norms (Fig 6)
    run('re100-cylinder/reconstruction_loop.m')
    
    %% Mixing layer
    % Reconstruction from downsampling (Fig 7)
    run('mixing-layer/ml_downsample.m')
    run('mixing-layer/ml_down_sparsity.m')
    
    % Reconstruction from windows (Fig 8)
    run('mixing-layer/ml_centerline.m')
    run('mixing-layer/ml_sparsity.m')
    
    %% NOAA sea surface temperature
    run('noaa-sst/remote_read.m')  % Download data from NOAA THREDDS server
    run('noaa-sst/sst_reconstruct.m')  % Reconstruction (Fig 9)
    
    %% HYCOM Gulf of Mexico surface vorticity
    run('hycom/remote_read.m');   % Download data from HYCOM THREDDS server
    run('hycom/hycom_reconstruct.m')  % Reconstruction (Fig 10)
    
    run('hycom/hycom_downsample.m')  % Downsampled reconstruction (Fig 11)
    
    run('hycom/kernel_demo.m')  % Demonstrate kernel scheme (Fig 14)
    
    %% Discussion (all data sets)
    % Error versus number of random point measurements (Fig 13)
    run('err_vs_sensors.m')  

end


%% Re=100 cylinder
SAVE_FIGS = false;

disp('Generating figure 1...')
GEN_FIG1;  % Components for diagram of reconstruction method

% Figure 2 is the flow chart (no code to generate)
% Figure 3 is four panels of the flow fields (also no code)

disp('Generating figure 4...')
GEN_FIG4;   % Singular value spectra for the different flows

disp('Generating figure 5...')
GEN_FIG5;  % Measurement strategies and err vs noise

disp('Generating figure 6...')
GEN_FIG6;  % err vs noise for different dictionaries and norms

pause
close all

%% Mixing Layer
% Reconstruction from downsampling
disp('Generating figure 7...')
GEN_FIG7;

% Reconstruction from centerline measurements
disp('Generating figure 8...')
GEN_FIG8;

pause
close all

%% NOAA SST
disp('Generating figure 9...')
GEN_FIG9;   % Example reconstructions of sea surface temperature

%% HYCOM
disp('Generating figure 10...')
GEN_FIG10;   % Example reconstructions of Gulf of Mexico vorticity

disp('Generating figure 11...')
GEN_FIG11;   % Example Gulf reconstructions from downsampling


pause
close all


%% Discussion figures (all data sets)
disp('Generating figure 12...')
GEN_FIG12;  % Uses output from pod_analysis.m (run above for figure 4)

disp('Generating figure 13...')
GEN_FIG13;  % Error vs number of random point measurements

%% Appendix: Generating kernels for HYCOM data
disp('Generating figure 14...')
GEN_FIG14;

disp('Done')