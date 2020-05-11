% The file 'fuelEconomy2.txt' contains fuel economy data for different models of cars (you can use the command edit fuelEconomy2.txt to examine its contents).

% In this exercise, you will import a subset of data and group the MPG (miles per gallon) values based on the number of cylinders and find the average MPG for each group. Finally, you will create a bar plot and customize it to look like the plot shown below:

% Step 1 – In this exercise, let's say we don't want to import the data for minivans and SUVs. You can do this by creating a datastore and importing only 362 lines from the file. Store the imported data into data.
 
% Step 2 – Group the variable NumCyl and then use the function splitapply to find the average CombinedMPG for each group. Store the result in avgMPG.
 
% Step 3 – Create a bar plot and customize it so that it looks like the plot shown above. Feel free to look for axes, bar, and figure properties in MATLAB Documentation. The figure and axes color should be set to [0.81 0.87 0.9] and the color of the bars should be set to [0,0.31,0.42]. For more help, you can click ‘Hint’ to get a list of properties that you will need to modify.

% Read data
dat = datastore('fuelEconomy2.txt');
dat.ReadSize = 362;
data = dat.read;

% Group by number of cylinders
[gNum,gVal] = findgroups(data.NumCyl);

% Find average by groups
avgMPG = splitapply(@mean,data.CombinedMPG,gNum);

% Create a bar chart
b = bar(avgMPG);
xlabel('Number of cylinders')
title('Average MPG')

% Customize the chart
f = gcf;
a = gca;

f.Color = [0.81 0.87 0.9];

a.Color = [0.81 0.87 0.9];
a.Box = 'off';
a.YAxisLocation = 'right';
a.YGrid = 'on';
a.GridColor = [1 1 1];
a.GridAlpha = 1;
a.XTickLabel = gVal;
a.YLim = [0 40];

ax = a.XAxis;
ax.TickDirection = 'out';

b.FaceColor = [0,0.31,0.42];
b.BarWidth = 0.5;





