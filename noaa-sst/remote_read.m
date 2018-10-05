% Download NOAA sea surface temperature data 1981-present from OpenDAP database
% Data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at http://www.esrl.noaa.gov/psd/
% code: Jared Callaham (2018)
%
% Matlab NetCDF/OpenDAP docs:
%  https://www.mathworks.com/help/matlab/import_export/importing-network-common-data-form-netcdf-files-and-opendap-data.html
clear all; close all; clc;

%% OpenDAP remote datasets
source_80s = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/noaa.oisst.v2/sst.wkmean.1981-1989.nc?time[0:426],sst[0:426][0:179][0:359]';
source_90pres = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/noaa.oisst.v2/sst.wkmean.1990-present.nc?time[0:1486],sst[0:1486][0:179][0:359]';
source_mask = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/noaa.oisst.v2/lsmask.nc?lat[0:1:179],lon[0:1:359],mask[0][0:179][0:359]';

lon = ncread(source_mask, 'lon');  % 'x' coordinates on map
lat = ncread(source_mask, 'lat');   % 'y' coordinates on map
mask = ~ncread(source_mask, 'mask');  % True for points on land
time_80s = ncread(source_80s, 'time');     % 1981-9
time_90s = ncread(source_90pres, 'time');  % 1990-present

%% Download data
sst = zeros(length(lat), length(lon), length(time_80s) + length(time_90s));

disp('Downloading 1981-1989')
idx = 1;
for i=1:length(time_80s)
    fprintf('%d/%d\n', i, length(time_80s));
    T = ncread(source_80s, 'sst', [1, 1, i], [Inf, Inf, 1]);  % Get temperature
    T(mask) = NaN;  % Mask land values
    sst(:, :, idx) = fliplr(T)';  % Flip upright for plotting    
    idx = idx+1;
end

disp('Downloading 1990-present')
for i=1:length(time_90s)
    fprintf('%d/%d\n', i, length(time_90s));
    T = ncread(source_90pres, 'sst', [1, 1, i], [Inf, Inf, 1]);  % Get temperature
    T(mask) = NaN;  % Mask land values
    sst(:, :, idx) = fliplr(T)';  % Flip upright for plotting 
    idx = idx+1;
end

sst = single(reshape(sst, [], size(sst, 3)));

%% Plot sample frame
% T = sst(:, end);
% 
% figure()
% pcolor(reshape(T, length(lat), length(lon)));
% shading interp; colormap(jet);


%% Save data
time = [time_80s; time_90s];
save('../data/sst_weekly.mat', 'lat', 'lon', 'sst', 'time', '-v7.3');

