%Provider	Connection Function
%Bloomberg®	% blp
%Thomson Reuters™ Datastream® API	% datastream
%FactSet®	% factset
%FRED®	% fred
%Haver Analytics®	% haver
%IQFEED®	% iqf
%Kx Systems®, Inc.	% kx
%SIX Financial Information	% tlkrs

% Connect to the server using a function from the table above;
c = fred; % Connect to FRED server;

% the data series required (e.g. 'USD12MD156N' to represent the daily 12-month London Interbank Offered Rate (LIBOR)),
% the start and end dates (specified as datetime variables).

tStart = datetime(2014, 1, 1);
tEnd = datetime(2014, 12, 31);
LIBOR = fetch(c, 'USD12MD156N', tStart, tEnd); % LIBOR
SP500 = fetch(c,'SP500', tStart, tEnd); % extract S&P500 index

% get usefulData from d, the Data part
Data_LIBOR = LIBOR.Data;
Data_SP500 = SP500.Data;

% Extract the first column of the matrix_dates
SP500_Dates = Data_SP500(:, 1);
SP500_Dates = datetime(SP500_Dates,'ConvertFrom', 'datenum'); % Convert these numeric dates to datetimes.

LIBOR_Dates = Data_LIBOR(:, 1);
LIBOR_Dates = datetime(LIBOR_Dates,'ConvertFrom', 'datenum'); % Convert these numeric dates to datetimes.

% Close the server connection and tidy up.
close(c)
clear c 

%% plot the Data
plot(LIBOR_Dates,Data_LIBOR(:,2));

hold on;
plot(Data_SP500(:,2));
hold off;

%% plot from common data

[commonDates,idxDAX,idxSP500] = intersect(daxDates,sp500Dates);
commonDAX = daxData(idxDAX,2);

commonSP500 = sp500Data(idxSP500,2);
commonData = [commonDAX commonSP500];

plot(commonDates,commonData)
plot(t, y, 'DatetimeTickFormat', 'QQQ MMM dd yyyy') % QQQ means Q1/2/3/4

plotyy(t1, y1, t2, y2) % Use plotyy to create a chart with two vertical axes

legend(['FTSE';'BARC'],'location','best','orientation','horizontal') % location and orientation of the legend



