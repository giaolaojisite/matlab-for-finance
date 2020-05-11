%% 1. spotCurve.m
% Fit data to spot curve data

%% Import Data
% Data is imported into a table named spotRates with two variables:
% maturity and spotrate
spotRates = readtable('spotRates.csv');

%% Plot the spot curve
figure
plot(spotRates.maturity,spotRates.spotrate,'o')
xlabel('Loan maturity (months)')
ylabel('Spot rate')
hold on

%% TODO: Fit a quadratic polynomial to the data
c = fitlm(spotRates,'quadratic');

%% TODO: Add the curve to the plot of the data
trendLine = predict(c,(1:15)');
plot(1:15,trendLine,'r')

%% TODO: Predict the rate for a maturity of 15-months called rate15
rate15 = trendLine(end);
title(['Quadratic model 15 month rate is ',num2str(rate15)])




%% 2. djiRegress.m
% Dow Jones Industrial Average Regression 
% using Macroeconomic Factors
%
% Practice performing regression to perform a
% fit based on multiple economic factor data to 
% predict DJI movement.

%% Import the data
djData = readtable('djiEconModel.csv');

%% TODO: Create a model to fit the data
model = fitlm(djData{:,3:end},djData.DJI);

%% TODO: Plot the original data with the fitted values in red
plot(djData.Dates,djData.DJI)
hold on
plot(djData.Dates,model.Fitted,'r','DatetimeTickFormat','QQQ-yy')
hold off
grid on
axis tight
title('Dow Jones Industrial Average')
legend('Original Data', 'Multiple Linear Fit', 'Location', 'NW')

%% TODO: Create Rsquared containing the ordinary R-squared value
Rsquared = model.Rsquared.Ordinary;
disp(['The Ordinary R-squared value is: ',num2str(Rsquared)])





%% 3. stocksDe.m - financial multivariate stats exercise

%% Import table of returns
returns = readtable('germanStocks.csv');

%% Plot original data
subplot(2,1,1)
plot(returns.Dates(end-50:end),returns.DAX(end-50:end))
hold on
grid on
title({'Deutscher Aktienindex (DAX) Returns';'German Stock Index Return Values'})

%% Set-up dimensions for Monte-Carlo Simulations
nSteps = 20;            % Number of steps into the future
nExperiments = 2;       % Number of different experiments to run

%% TODO: Generate random numbers from the information in DAX
tFit = fitdist(returns.DAX,'tLocationScale');
simReturns = random(tFit,nSteps,nExperiments);

%% Plot the returns
futureDates = returns.Dates(end) + days(1:nSteps);
plot([returns.Dates(end) futureDates],[returns.DAX(end)*ones(1,nExperiments) ; simReturns])
hold off

%% Plot the Index Values
lastDAX = 6147.97;
predictions = ret2tick(simReturns,'StartPrice',lastDAX);
subplot(2,1,2)
plot([returns.Dates(end) futureDates],predictions)
grid on
title({'Deutscher Aktienindex (DAX)';'German Stock Index'})

