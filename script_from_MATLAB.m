% Script to predict future prices using DBA model

% Load data and identify the stock of interest
load stockData
stock = BA;

% Calculate the daily log returns and the statistics
stockReturns = diff(log(stock));
mu = mean(stockReturns);
sigma = std(stockReturns);

% Simulate prices for 22 days in future based on GBM
deltaT = 1;
S0 = stock(end);
epsilon = randn(22,200);
factors = exp((mu-sigma^2/2)*deltaT + sigma*epsilon*sqrt(deltaT));
lastPriceVector = ones(1,200)*S0;
factors2 = [lastPriceVector;factors];
paths = cumprod(factors2);

% Extract final price from paths
finalPrices=paths(end,:);

% Log returns 
possibleReturns=log(finalPrices/S0);

% plot histogram
histogram(possibleReturns,20);
ylim([0,35]);
% line for the 5 percentile and plot
var5=prctile(possibleReturns,5);
hold on;
plot([var5,var5],[0,35],'-r','LineWidth',1);
hold off;

%% answer from MATLAB
% Load data and identify the stock of interest
load stockData
stock = BA;

% Calculate the daily log returns and the statistics
stockReturns = diff(log(stock));
mu = mean(stockReturns);
sigma = std(stockReturns);

% Simulate prices for 22 days in future based on GBM
deltaT = 1;
S0 = stock(end);
epsilon = randn(22,200);
factors = exp((mu-sigma^2/2)*deltaT + sigma*epsilon*sqrt(deltaT));
lastPriceVector = ones(1,200)*S0;
factors2 = [lastPriceVector;factors];
paths = cumprod(factors2);

% Extract final prices
finalPrices = paths(end,:);

% Calculate returns
possibleReturns = log(finalPrices) - log(S0);

% Plot the histogram of returns
histogram(possibleReturns,20)

% Calculate the VaR
var5 = prctile(possibleReturns,5);
hold on
plot([var5 var5],[0 40],'r')
hold off
