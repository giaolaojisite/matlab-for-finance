%% TODO - Read in data from stockData.csv and stockInfo.csv
% Create a file named stocks.csv with the stock name and mean value

%% Read in data
stockData = readtable('stockData.csv'); % date in a column and each column represents a company
stockInfo = readtable('stockInfo.csv'); % company names in a column

%% Find average prices
avgRet = varfun(@mean,stockData(:,2:end)); % exclude the first column of date

%% Create a table with just company names and average prices
stocks = stockInfo(:,1);
stocks.MeanPrices = avgRet{1,:}';

%% Write the data out to a file named stocks.csv
writetable(stocks,'stocks.csv')
