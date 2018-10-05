% Download HYCOM Gulf of Mexico data 1992-present from THREDDS database
% source: https://hycom.org/
% code: Jared Callaham (2018)
%
% Matlab NetCDF/OpenDAP docs:
%  https://www.mathworks.com/help/matlab/import_export/importing-network-common-data-form-netcdf-files-and-opendap-data.html
clear all; close all; clc;
 
%% Load expt_19.0 (10/2/92-7/31/95)
source_92 = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0?lat[0:1:2000],lon[0:1:4499],time[0:1:1030]';

lon = ncread(source_92, 'lon');  % 'x' coordinates on map
lat = ncread(source_92, 'lat');   % 'y' coordinates on map
t_92 = ncread(source_92, 'time');

lon_idx = find( (lon >= -98) & (lon <= -76.4) );
lat_idx = find( (lat >= 18.1) & (lat <= 32));
lat = lat(lat_idx); lon = lon(lon_idx);

% % Redefine source to include velocities
source_92 = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0?water_u[0:1:1030][1][1227:1:1400][1025:1:1295],water_v[0:1:1030][1][1227:1:1400][1025:1:1295]';
 
idx = 0;
for j=1:length(t_92)
    idx = idx + 1;
    fprintf('Loading %d/%d\n', j, length(t_92))
    uvel(:, :, idx) = (ncread(source_92, 'water_u', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % East-west velocity
    vvel(:, :, idx) = (ncread(source_92, 'water_v', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % North-south velocity
end

%% Load expt_19.1 (8/1/95-12/31/2012)
source_95 = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1?time[0:1:6324],water_u[0:1:6324][1][1227:1:1400][1025:1:1295],water_v[0:1:6324][1][1227:1:1400][1025:1:1295]';
t_95 = ncread(source_95, 'time');

for j=1:length(t_95)
    idx = idx + 1;
    fprintf('Loading %d/%d\n', j, length(t_95))
    uvel(:, :, idx) = (ncread(source_95, 'water_u', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % East-west velocity
    vvel(:, :, idx) = (ncread(source_95, 'water_v', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % North-south velocity
end

%% Load expt_90.9 (5/12-8/20/13)
% Note: start later to throw away overlapping data (and lon changes)
source_12 = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_90.9?time[259:1:460],water_u[259:1:460][1][1227:1:1400][3275:1:3545],water_v[259:1:460][1][1227:1:1400][3275:1:3545]';
t_12 = ncread(source_12, 'time');

for j=1:length(t_12)
    idx = idx + 1;
    fprintf('Loading %d/%d\n', j, length(t_12))
    uvel(:, :, idx) = (ncread(source_12, 'water_u', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % East-west velocity
    vvel(:, :, idx) = (ncread(source_12, 'water_v', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % North-south velocity
end

%% Load expt_91.0 (8/21/2013 - 4/4/2014)
source_13 = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.0?time[0:1:233],water_u[0:1:233][1][1227:1:1400][3275:1:3545],water_v[0:1:233][1][1227:1:1400][3275:1:3545]';
t_13 = ncread(source_13, 'time');

for j=1:length(t_13)
    idx = idx + 1;
    fprintf('Loading %d/%d\n', j, length(t_13))
    uvel(:, :, idx) = (ncread(source_13, 'water_u', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % East-west velocity
    vvel(:, :, idx) = (ncread(source_13, 'water_v', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % North-south velocity
end

%% Load expt_91.1 (4/4/2014-4/18/2016)
source_14 = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1?time[0:1:737],water_u[0:1:737][1][1227:1:1400][3275:1:3545],water_v[0:1:737][1][1227:1:1400][3275:1:3545]';
t_14 = ncread(source_14, 'time');

for j=1:length(t_14)
    idx = idx + 1;
    fprintf('Loading %d/%d\n', j, length(t_14))
    uvel(:, :, idx) = (ncread(source_14, 'water_u', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % East-west velocity
    vvel(:, :, idx) = (ncread(source_14, 'water_v', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % North-south velocity
end

%% Load expt_91.2 (4/18/2016-present)
source_16 = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2?time[0:1:737],water_u[0:1:737][1][1227:1:1400][3275:1:3545],water_v[0:1:737][1][1227:1:1400][3275:1:3545]';
t_16 = ncread(source_16, 'time');

for j=1:length(t_16)
    idx = idx + 1;
    fprintf('Loading %d/%d\n', j, length(t_16))
    uvel(:, :, idx) = (ncread(source_16, 'water_u', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % East-west velocity
    vvel(:, :, idx) = (ncread(source_16, 'water_v', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % North-south velocity
end

%%
uvel = single(uvel); vvel=single(vvel);
t = [t_92; t_95; t_12; t_13; t_14; t_16];
%save('../data/gom_00Z.mat', 'uvel', 'vvel', 'lat', 'lon', 't', '-v7.3');  % File is ~1.3 Gb

[ny, nx, m] = size(uvel);
n = nx*ny;

%% Convert to vorticity and save
vort = zeros(n, m);
for i=1:m
    [omega_z, ~]= curl(X,Y,uvel(:, :, i),vvel(:, :, i));
    vort(:, i) = omega_z(:, :, 1);
end

clearvars -except vort lat lon t
save('../data/gom_vort.mat')  % File will be ~1.1 Gb