%% Import data
data = readtable('marchTemps.csv');
%% Extract data
lon = data.longitude;
lat = data.latitude;
T = data.temperature;
% Create vectors for new longitude and latitude values
lonvec = 110:155;
latvec = -45:-10;
% Create a grid of longitudes and latitudes
[longrid,latgrid] = meshgrid(lonvec,latvec);
% Interpolate the scattered March temperature data onto the grid
Tgrid = griddata(lon,lat,T,longrid,latgrid);
%% TODO: Create two images of the temperature data
figure
im = pcolor(longrid,latgrid,Tgrid);
im.EdgeColor = 0.5*[1 1 1];
im.EdgeAlpha = 0.5;

figure
im2 = imagesc(lonvec,latvec,Tgrid);
axis xy
