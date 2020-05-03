% The function of the script

help filename % return the function/description of the script

doc sum %% search function 'sum'
help sum %% show directly

edit filename.m %% create a new script/open the script
filename %% run the script

clear variable_name; % clear the specific variable in the workspace
disp(tablename) % display the table

save('filename.mat','x','y');
load('filename.mat');
plot
class(x)  %% show the type of x, like numeric, character,

docsearch('set operations') % search in MATLAB
openExample('finance/ConvertPriceSeriesToReturnSeriesTimeTableExample')

plot %% default color-blue default plot-line

dataset=(x1:interval:x2); with' means column vector;
or x=[x1:x2] %%better or x=[x1,x2,x3,x4], or x = x1:x2,

assetNames = {'BARC' 'GSK' 'BRBY' 'SVT'}

dataset=linspace(x1,x2,numbers between); %% default with 100 numbers;

BARCsubset=BARC([1:9]); 

[minFTSE,idx] = min(FTSE); %% find the location of max/min number;

var=prctile(BARCReturns,5); %% find the x interception 'var' at 5%

%% plot 
subplot(2,1,1) %% subplot with a plot,
subplot(4,4,[3,8]);
subplot(4,4,[3,4],'clear');

plot(FTSE,'LineWidth',3,'Color',[0.3 0.3 0.8]) %% between 0 and 1 representing the red, green, and blue components of the color desired.
xlim([1 100]); 

'MarkerSize' default 6, 
'MarkerFaceColor', 
'MarkerEdgeColor',
'LineWidth', default 0.5,

xlable('xtitle'),  %% '\beta', '\sigma';
ylable('ytitle'),
title('title'),      
title(['Value at risk = ' num2str(Value)]);
title(['title in ' char(number)]), %% char(163)=£, char(36)=$; use[] to combine string into title,

legend('y=x','location','best'), %% southwest, north etc.,

histogram(dataset,'Display','stairs',),
histogram(GSKRetns,40,'EdgeColor','none')
plot(data,number of bars,'LineWidth',3,'MarkerEdgeColor','r')

plot([x1,x2],[y1,y2]) %%
figure(1), %% creat a figure without closing the present figure,

hold on, %% plot on the same figure,
hold off, 

%% mathworks.com/help/matlab/ref/plot.html%namevaluepairarguments
%% https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.histogram.html%namevaluepairarguments 

[Principal,Interest,Balance,Payment] = amortize(Rate,NumPeriods,PresentValue);

output = M([2 3],[1 2 3 4]); or 1:4; or [1:4]; or 1:end, or : if it is the entire column/row;
output = M(end-1,[2]);
output = M([end-2,end-1],[end-1,end]);

[num,txt] = xlsread('stockData','Prices')

round(dataset,number of digits); 
round(dataset,number,'significant');

mean(dataset,1) or mean(dataset) %%% the results for each column;
diff(dataset,1) or diff(dataset)
min(dataset)
std(dataset) 

mean(dataset,2) %%% the results for each row;
min(dataset,[],2)
std(dataset,??) 


factors=exp((mu-sigma^2/2)*deltaT+sigma*epsilon*sqrt(deltaT));

%% concatenate matrix
m3 = [m1 m2] or [m1,m2] %% horizontal concatenating
m3 = [m1;m2] %% vertically concatenating

B = cumprod(A) %% default horizontally, 2 means vertically cumproduction;

%% date 
d = datetime(1999, 12, 31); % t = 21-Oct-2005
t = datetime(2005,10,21,17,10,0); % t = 21-Oct-2005 17:10:00
t = datetime('Oct 21, 2005');
t = datetime('21 Oct 2005','InputFormat','dd MMM yyyy'); % MMM represents Oct
d=datetime('2015-09-30','InputFormat','yyyy-MM-dd'); % MM represents 09
>> d = datetime('July 14 2015', 'InputFormat', 'MMMM dd yyyy')
>> s = {'Q1 2015', 'Q2 2015'}
>> d = datetime(s, 'InputFormat', 'QQQ yyyy') % QQQ
d = 01-Jan-2015   01-Apr-2015

d = datetime(1990, (1:12).', 1) % or 
d = datetime(1990, (1:12), 1).' % or
d = datetime(1990, (1:12), 1)' 
d = linspace(datetime(1995,1,1),datetime(1995,6,30),20); % find the dates in between
d = d1:interval_in_days:d2; % find the dates in between
w = (t1:calquarters(1):t2).' % only in this way for calendar duration

t = size(datetime(2010, 8, 1):datetime(2010, 10, 1)); % measure the dimenstion

d.Format % check the format
d.Format = 'eee MM yyyy G' % change the format, e-Fri, G-CE;
d = datetime(2008,11,22,'Format','eee MMM dd yy') % directly change the format

Format=d.Format; % store the format;
d = datetime(2008,11,22,'Format',Format); % use the stored format;

t = duration((1:5).', 0, 0) % default format h m s;
% quick means of duration;
t = years(1:5); % a 1 by 5 array;
t = minutes(0:10:60) % with space of 10mins;

t = calendarDuration(1,2,3,4,5,6) % year, month, and day, plus the hour, minute, and second;
% quick means of calendar duration;
t = calyears(2) % 2 calendar years;

LIBOR_Dates = datetime(LIBOR_Dates,'ConvertFrom', 'datenum') % Convert these numeric dates to datetimes.
% convert from date numbers to date time

[c,idxA,idxB] = intersect(A,B) % idx represent the location of intersection in A/B;
commomDAX=daxData(c,2); % find the correspondant value at the common dates;

plot(t, y, 'DatetimeTickFormat', 'QQQ MMM dd yyyy') % or, QQQ means Q1/2/3/4
plot(t, y, 'DurationTickFormat', 'MMM dd yyyy')

dy = day(d, 'dayofyear') % extract day of the year
dy = year(d); % d-datetime, 'year' function will read the year of d; 'years' function measurea duration by num of years
t = datetime(2005, 10, 21);
[yr,mth,d] = ymd(t) % extract the year etc. from data
[hr,mn,sd] = hms(t) % extract the hour etc. from data

cqr = split(c, 'quarters') % 'years' extract from calendarDuration
[cmth, cwk] = split(c, {'months', 'weeks'})

d=calweeks(2); % convert number to duration;
d=calweeks(d); % convert duration back to number (how many weeks does it include);
v=t1+days(1:15).'; % a column vector

dt=datetime1-datetime2 % the interval in terms of duration, days, hours, seconds;
nMths = numel(t1:calmonths(1):t2) - 1 % datetime substraction=duration, not calendar duration,
nDays=dt/days(1); % the number of days in between

cw = caldiff(d1, {'weeks', 'days'}) % the time difference in weeks/days
cd = between(d2, t2) % format year, month, day
d1Mon = dateshift(d1, 'dayofweek', 'Monday', 'nearest') % 0 means nearest, 1 means current or next, -1 means current or before, 2 means after current or next
d2MonthEnd = dateshift(d2, 'end', 'month', 'nearest') % end/start of month/week
t = datetime('today') 
t2 = dateshift(t,'start','month','next')


%% import to table

stockPrices = array2table(Prices,'VariableNames',{'BA','BARC','BRBY'}); % convert array/matrix to table with titles,
stockData = readtable('stockPrices.xls') % read data from excel
% .txt, .dat, or .csv for delimited text files
% .xls, .xlsb, .xlsm, .xlsx, .xltm, .xltx, or .ods for spreadsheet files
brby = stockData(5:9,[1,4]) % extract part from the table
BA = stockData(:,'name_of _column') % extract from the column 'BA' with ()
rStocks = stockData(:,{'RBS','RDSA'}) % extract from the columns with ()

stockData.BRBY % the numbers of the column
minNovBrby = min(stockData.BRBY(7:end)) % filename.variablename to extract the number of that column
>> baRbs1 = stockData{:,[2 5]}; % extract number with '{}'
>> baRbs2 = stockData{:,{'BA','RBS'}}; % extract number with '{}'

stockData.BARC=709*ones(height(stockData.BARC),1) % change number of a column to 709
height(table) % measure the num.rows in a table
>> stockData.FTSE = ftse % add a column (from data ftse) with the name FTSE, 
>> stockData.RBS = round(stockData.RBS) % change the column
>> stockData.BARC = 709*ones(height(stockData),1)

stockData.BA=stockData.BA.*stockData.ExchangeRate; % change the number of single column;
stockData{:,2:end} = 1.87*stockData{:,2:end}; % modify data in Table with '{}';
stockData{:,2:end} = round(stockData{:,2:end}) % round the numbers
stockData.Dates = datenum(stockData.Dates) % change date type to date numbers;

stockProps = stockData.Properties % save the table properties;
names = stockData.Properties.VariableNames % save a specific property 
stockData.Properties.VariableNames{1}='DateNum'; % change the variable name with '{}', or '{'Date'}';

% function handle
f = @mean;
fx=f(x); % mean of x;

x = fzero(@fcnName,xStart) % find x when y=0,xstart is an initial start point for the function to find;
x = fminbnd(@fcnName,lowerBnd,upperBnd) % find x when y approaches min within the x range

y = varfun(@afctname,x) % apply function to the table x, results stored in table y

% categorical
stockData.Industry=categorical(stockData.Industry) % organize repeated text of a variable
x = categorical([ 2 2 1 2 3 1 ]);
x = mergecats(x,{'2','3'},'C');

upcycle = {'Consumer Cyclical','Services','Transportation','Financial'};
downcycle = {'Healthcare','Capital Goods','Energy','Utilities','Basic Materials'};
stockInfo.Industry = categorical(stockInfo.Industry) % categorical the column
stockInfo.Classification = stockInfo.Industry % create a new column
% change the name of some to upcycle
stockInfo.Classification = mergecats(stockInfo.Classification,upcycle,'upcycle') 
stockInfo.Classification = mergecats(stockInfo.Classification,downcycle,'downcycle')
boxplot(stockInfo.MeanReturns,stockInfo.Classification) % boxplot the returns by categories,
% create a box plot of data in X that is grouped by G,
boxplot(X,G);

% export tables
writetable(tableName,'filename.xls','Range','C5:L23') % Specify the spreadsheet range.
writetable(tableName,'filename.txt','Delimiter','\t') % Specify the text delimiter.
writetable(stockInfo,'stocksTab.txt','Delimiter','\t'); % tab-delimited

% trading stragegy

expMovAvg = movavg(FTSE,'exponential',n); % exponential moving average of n samples
%%% moving average convergence/divergence (MACD)
%%% leading (short-term) exponential moving average -  lagging (long-term) exponential moving average
stockPrices = readtable('stockData.csv');
FTSE = stockPrices.FTSE;
% TODO - Calculate the leading and lagging moving averages and then calculate the MACD.
movAvgShort = movavg(FTSE,'exponential',3);
movAvgLong = movavg(FTSE,'exponential',5);
MACD = movAvgShort - movAvgLong;

% logical operator, >, <, >=, <=, ==, and ~=
% '&'and,'|'or,'~'not,
test = ~(1 < 2); % not the result of 1<2
test = pi>3&pi<3.2; 

I = v<0.005; % check an array and store the results in an array;
check1 = any(I); % check whether any of the elements of a logical vector are true
check2 = all(I); % check whether all of the elements of a logical vector are true
countSmall=nnz(I); % nnz (number of non-zeros) to compute the number of elements of an array that are not zero
test1 = 'abc' == 'abd'; % results is logical array 110
test2 = strcmp('abd','abc'); % string comparison, result is logical 0
test3 = strcmp(tickers,'BARC'); % compare a string with multiple strings stored in a cell array 

%%% MACD
stockPrices = readtable('stockData.csv');
FTSE = stockPrices.FTSE;
movAvgShort = movavg(FTSE,'exponential',3);
movAvgLong = movavg(FTSE,'exponential',5);
MACD = movAvgShort - movAvgLong;
% num of days in 2007 with MACD>=5
up=nnz(MACD>5);
stockPrices.MACD=MACD;
up2007=nnz(year(stockPrices.Dates)==2007&stockPrices.MACD>=5);
% TODO: Create a column vector named econPerformance whose elements
% are 1 when the economy is up, -1 when the economy is down, and 0 otherwise.
econPerformance = zeros(length(MACD),1); % or smarter
econPerformance = 0*MACD; % the same dimension as MACD
econPerformance(MACD >= 5) = 1;
econPerformance(MACD <= -5) = -1;
% Plot 
plot(Dates,econPerformance)
ylim([-2,2])

%%% trading strategy project
% Load data
stockPrices = readtable('stockData.csv');
stockInfo = readtable('stockInfo.csv');
% Determine economy's Performance
movAvgShort = movavg(stockPrices.FTSE,'exponential',3);
movAvgLong = movavg(stockPrices.FTSE,'exponential',5);
MACD = movAvgShort - movAvgLong;
econPerformance = zeros(length(MACD),1);
econPerformance(MACD >= 5) = 1;
econPerformance(MACD <= -5) = -1;
stockPrices.econPerformance = econPerformance;
% Examine the variable stockPrices - the momentum on Nov-6-1986 is 1 i.e., 'up'.
% Extract the stock codes of all the 'upcycle' stocks.
upcycleIdx = strcmp(stockInfo.Classification,'upcycle');
stocksToBuy = stockInfo{upcycleIdx,'Code'};
% Extract the corresponding prices
stocksToBuyPrices = stockPrices{:,stocksToBuy};

% extract data using logical indexing
I = v<0.005;
r = v(I); or
r = v(v < 0.005); % in one step
v(v < 0.005) = 0; % replace values with logical indexing

% deal with missing values
xAvg = mean(x,'omitnan'); % ignore missing value, works for functions 'mean','std','median','cov','var'
% 'max','min' by default ignore N/A values
sdmat = stockData{:,:}; % store the numerical data
stockReturns = tick2ret(sdmat); % calculate the return (like diff)

%%% find and display
%TODO: find minHL, maxHL, and avgHL
minHL = min(stocks.HL);
maxHL = max(stocks.HL);
avgHL = mean(stocks.HL,'omitnan');
% display 
disp(['Minimum Value of HL: ',num2str(minHL)]) % Minimum Value of HL: 177
disp(['Maximum Value of HL: ',num2str(maxHL)])
disp(['Average Value of HL: ',num2str(avgHL)])

test = isequal(x,y); % determine if x=y, will return 0 for table containing NaN
test = isequaln(x,y); % determine if x=y, will return 1 for table containing NaN

% locating missing data
x = [ 1 Inf 13 NaN pi ]; 
isnan(x); % determine if NaN, NaN with 1, and all the other with 0, 'isnan' - 'ismissing'
isfinite(x); % determine if finite
x(isnan(x)) = 1; % replace the NaN value with 1
x(isnan(x)) = 1; % delete NaN value (the array will get shorter)
% remove rows with missing data
>> nnp = nnz(isnan(data.powers));
>> nnd = nnz(isnat(data.dates));
>> isHole = ismissing(data); % determine the whole dataset with all types
>> holeRow=any(isHole,2); % the row with any missing data, 2 means scaning row by row
>> fullData=data(~holeRow,:) % only contain rows with data using '~' 

x = table((1:10)',rand(10,1)); % create a 10 by 1 table

% find the number of rows with missing data
>> totRowsNM=nnz(~(any(ismissing(stockData),2))) % or remove rows and use 'height' to measure the new row, 'any' - 'all'

% import data, the headings get stored as variable names(dont count rows)
>> stockData = readtable('stockDataNanCol.txt')
% calculate returns with 'tick2ret', 1 less row than the original table
>> stockRet = tick2ret(stockData{:,:})
% determine which stock loses all data
>> missingStocks = all(ismissing(stockData),1) % determine stock by stock (column by column)
% read the variable names from a table
>> stockCodes = stockData.Properties.VariableNames 
% read the variable name with missing data
>> emptyStocks = stockCodes{missingStocks}

x(:,all(ismissing(x))=[]; % remove the entire column

% this might cause trouble, think about the sequence of removing data
% think about how to preserve the most data
stockRet = stockRet(~any(ismissing(stockRet),2),~any(ismissing(stockRet)));

%%% remove all NaN
stockData = readtable('stockDataNanCol.txt');
% Calculate returns
stockRet = tick2ret(stockData{:,:});
% TODO: Remove NaNs from stockRet
stockRet(:,all(isnan(stockRet))) = [];
stockRet(any(isnan(stockRet),2),:) = [];
% Show the return values
disp('stockRet = ')
disp(stockRet)

% interpolating missing data
idx = ~isnan(y) % y with valid data
yClean = y(idx);
xClean = x(idx); % use X to predict missing Y
y_predict = interp1(xClean,yClean,x(~idx),'previous','extrap') % with 'extrap' to extrapolate data (predict outside the range)

% with x,y, want to plot xfine and yfine (denser plot)
xFine = 0:0.1:13;
yFine = interp1(x,y,xFine,'spline'); % seems 'spline' will auto extrapolate

%%% table interpolate with function
function y = tableInterp(y)
% Create a vector of evenly spaced x-values the same length as y
x = (1:length(y))';
% Create a logical vector with true at the location of non-missing values
idx = ~isnan(y);
% Extract non-missing values
xClean = x(idx);
yClean = y(idx);
% Use this index to extract non-missing data and interpolate the missing
% values
y = interp1(xClean,yClean,x,'spline','extrap');
end % end function
% interpolate table t
tInterp = varfun(@tableInterp,t);

% 'corrcoef' determine the matrix of pairwise correlation coefficients between each pair of columns
cRaw = corrcoef(pricesRaw);
% Visualize the data with imagesc
imagesc(cRaw);

%%% plot denser data (finer data points)
% Interpolate currency exchange rate data
% IMPORT DATA from CanadianExchange.csv
canEx = readtable('CanadianExchange.csv');
canEx.Dates = datetime(canEx.Dates);
% TODO: Plot the Canadian dollar exchange rates using black points
plot(canEx.Dates,canEx.ExchangeRate,'k.','MarkerSize',6)
% Interpolate the missing values
% Create a duration vector containing the current dates in the plot 
durationfine = min(canEx.Dates):max(canEx.Dates);
% Interpolate the missing values
exfine = interp1(canEx.Dates,canEx.ExchangeRate,durationfine);
% Plot the interpolated data in red
hold on
plot(durationfine,exfine,'r')
xlabel('Dates')
ylabel('Canadian Exchange Rate')
grid on

%%% make the dataset symmetrical
% TODO: Replace NaNs in indCorr to make it symmetric
indCorrT=indCorr';
indCorr(isnan(indCorr))= indCorrT(isnan(indCorr))

%%% use function to find largest N/smallest N of a dataset
% Do not edit %
x = rand(1,10);
%%%%%%%%%%%%%%%
% TODO: Modify the function call
top5 = getLargestN(x,5);
top7 = getLargestN(x,7);
bottom5 = getSmallestN(x,5);
% TODO: Modify the function definition
function t = getLargestN(v,N)
sortedV = sort(v,'descend');
t = sortedV(1:N);
end
function b = getSmallestN(v,N)
    sortedV=sort(v,'descend');
    b=sortedV(end-N+1:end);
end

%%% simulates the future stock prices using the historical prices provided in the input and outputs the N% VaR. 
function [var,finalPrices] = calculateVaR(stock,N)
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
% Calculate the VaR
var = prctile(possibleReturns,N);
end

% set path to run matlab code in other folders
var = calculateVaR(BARC,5);
varianceBARC = var(BARC); % will cause conflicts
% 'which' function to check the precedence
which var; % show precedence
which var -all; % show precedence and location

% mortgageCall function

% flow control
data=input('Enter number of hours: '); % expect a numeric input
name=input('Enter the employee name: ','s'); % expect a string input as mentioned by 's'

disp('MATLAB Rocks') % Display text or array
fprintf('Cost = $%d\n',20)	% Write data on the screen
warning('Range is out of bounds.')	% Display warning message in orange, program continues to execute
error('Data not found.') % Display error message in red, program stops execution

% if statement
if test
  statements 1
end
% if-else statement
if test
  statements 1
else
  statements 2
end
% if-elseif-else statement
if test1
  statements 1
elseif test2
  statements 2
elseif test3
  statements 3
else
  statements 4
end

% use ==,not = in the if test

%%% return stocks to buy with performance
% Load data
stockPrices = readtable('stockData.csv');
stockInfo = readtable('stockInfo.csv');
% Calculate the economy's performance
% for a specific date
dateOfInterest = datetime([2006 11 22]);
performanceAllDates = calculateEconPerformance(stockPrices.FTSE);
performance = performanceAllDates(stockPrices.Dates == dateOfInterest);
if performance == 1
   upcycleIdx = strcmp(stockInfo.Classification,'upcycle'); % compare the columns, find rows matching upcycle
   stocksToBuy = stockInfo{upcycleIdx,'Code'}; % the corespondant rows in the 'codes' column
elseif performance == -1
   downcycleIdx = strcmp(stockInfo.Classification,'downcycle');
   stocksToBuy = stockInfo{downcycleIdx,'Code'}; % mind the '{}'
end

% Run the switch-case code.	
x = 2;
switch x
    case 1
        disp('x is 1')
    case 2
        disp('x is 2')
    otherwise
        disp('x is neither 1 nor 2')
end
x is 2 % The output shows that the statements in the case where x==2 are executed.	

% while 
while x> 0.1 % excecute until x<=0.1
    disp(x)
    x = x-1;
end

% calculate the corelation between two columns of the array
c = corr(M)
c = 
    1.00   0.85
    0.85   1.00
real_c=c(1,2) % real_c=c(1,2)

%%% corelation and loop every 15 rows 
% using the function corr. The corr function outputs a matrix.
% Extract the correlation coefficient from the output matrix and store it in c.
cm = corr(indexRetns);
c = cm(1,2);
% Task 3 - Write a 'for loop' to compute the correlation of 15 elements at a time.
% Store the correlation coefficients in a vector called rollingCorr.
windowSize = 15;
numRecords = size(indexRetns,1);
for k = 1:(numRecords - windowSize + 1)
   rollingCorrMatrix = corr(indexRetns(k:k+windowSize-1,:)); 
   rollingCorr(k) = rollingCorrMatrix(1,2);
end

% fitting models to empirical data with 'fitlm' function
>> FTSEfit = fitlm(dates,FTSE); % will not work because dates is not numeric/categorical array
>> d = days(dates-dates(1));
>> FTSEfit = fitlm(d,FTSE); % default linear or
>> FTSEfit = fitlm(table); % combine these as a matrix, or
>> FTSEfit = fitlm(FTSEtable,'PredictorVars','days','ResponseVar','FTSE'); % specify the variables in the table
% access the prediction properties
>> c = FTSEfit.Coefficients
>> rs = FTSEfit.Rsquared %  coefficient of determination, R2, gives an interpretation of goodness-of-fit.
>> f3 = FTSEfit3.Formula % The Formula property stores the mathematical representation of the regression model
%
>> FTSEfit2 = fitlm(d,FTSE,'quadratic');
>> FTSEfit3 = fitlm(d,FTSE,'poly3'); % a 3-degree linear regression model
plot(FTSEfit2); % plot the fit data, or
plot(d,FTSEfit.Fitted,'g'); % Fitted property of the LinearModel object stores the predicted response values
% predict the model performance in the future
>> y = predict(modelName,(1:300)'); % for example the model runs from 1:200, then 'predict' function predict 1:300;
>> y = predict(FTSEfit,(0:200)'); % predict the FTSEfit model

% residual

>> FTSEfit = fitlm(daysNum,FTSE);
>> plot(dates,FTSE);
>> hold on
>> plot(dates,predict(FTSEfit,daysNum))
% A normal distribution is one indication of a good model
plotResiduals(FTSEfit) % ranges of the residuals and their frequencies.
% the residual plot should not contain any predictive information and you should see a more symmetric and constant spread throughout the range.
plotResiduals(FTSEfit,'probability') %how the distribution of the residuals compares to a normal distribution with matched variance.
plotResiduals(FTSEfit,'lagged') %shows the serial correlation among the residuals.
plotResiduals(FTSEfit,'fitted') %shows the size of residuals versus the size of the fitted data.
plotResiduals(FTSEfit,'symmetry') %shows the where the residuals are located with respect to the median value
plotResiduals(FTSEfit,'caseorder') %shows the size of residuals versus the row order.

% lillietest on the residual data (one column from the model)
>> resTable = FTSEfit.Residuals
>> res = resTable.Raw;
>> [hLil,pLil] = lillietest(res) % two outputs hLil,pLil
% Since the hLil value is 1, this implies that the test rejects the null hypothesis that the data comes from a normal distribution. 
% A small p-value, under 5%, also supports the rejection of the null hypothesis.

% The 'fitdist' function creates a MATLAB probability distribution object.
>> pd = fitdist(x,distname) % a type of distribution of x
>> ndo = fitdist(FTSEreturns,'Normal');
>> ndo = fitdist(FTSEreturns,'Normal');
% create a t Location-Scale probability distribution object named tFit
>> tFit = fitdist(FTSEreturns,'tLocationScale')
% icdf function to determine the inverse cdf value of the fit at a given value.
>> pv95 = icdf(tFit,0.05); % x value at the point (5%)

% assign the property 'Normalization' to a value of 'pdf' in order to normalize the height of each bar such that the sum of the bar areas is 1.
>> histogram(FTSEreturns,binEdges,'Normalization','pdf')

% use the pdf function to determine the pdf value of the fit with a given input vector.
>> r = -0.05:0.001:0.05;
>> p = pdf(tFit,r);
>> plot(r,p,'r');

% generating random numbers
>> randomVector = random(tFit,1,4); % tFit is a probability distribution, 1 by 4 array;

% 'rng' function
% Seed the random number generator using 13.
>> rng('shuffle') % will use the current time to set the seed, avoid repeating

% t Location-Scale probability distribution object

% fit and predict future










































































