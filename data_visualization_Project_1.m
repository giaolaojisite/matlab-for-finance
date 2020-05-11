% The file 'fuelEconomy.txt' contains fuel economy data for different models of cars (you can use the command edit fuelEconomy.txt to examine its contents).

% The goal of this exercise is to divide the combined miles per gallon data into three classes and create a scatter plot of city and highway miles per gallon like the one shown below.

% Follow the steps listed below:

% Step 1 – Import the data into a table.

% Step 2 – The variable CombinedMPG contains missing values represented by NaN. Remove the rows corresponding to these missing values from the imported table.

% Step 3 – Discretize the CombinedMPG variable into three classes called 'Low', 'Medium', and 'High' using the bin edges 0, 20, 30, and 70 (values greater than or equal to 0 and less than 20 are classified as 'Low' and so on). Store the discretized result in a variable named MPGClass.
 
% Step 4 – Create a scatter plot with CityMPG on the x-axis and HighwayMPG on the y-axis. Observations with 'Low' combined MPG should be colored red ('r'), 'Medium' should be blue ('b'), and 'High' should be 'black' ('k').

% Read data
data = readtable('fuelEconomy.txt');

% Identify the rows containing NaN values and remove them
nanIdx = ismissing(data.CombinedMPG);
data(nanIdx,:) = [];

% Discretize Combined MPG
MPGClass = discretize(data.CombinedMPG,[0 20 30 70],{'Low' 'Medium' 'High'});

% Optional - Convert the class into a categorical array
MPGClass = categorical(MPGClass);

% Extract observations for various classes and plot them
scatter(data.CityMPG(MPGClass == 'Low'),data.HighwayMPG(MPGClass == 'Low'),...
'r','filled')
hold on
scatter(data.CityMPG(MPGClass == 'Medium'),data.HighwayMPG(MPGClass == 'Medium'),...
'b','filled')
hold on
scatter(data.CityMPG(MPGClass == 'High'),data.HighwayMPG(MPGClass == 'High'),...
'k','filled')
hold off

grid on
xlabel('City MPG')
ylabel('Highway MPG')
legend('Low Combined MPG','Medium Combined MPG','High Combined MPG');