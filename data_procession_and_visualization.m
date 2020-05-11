% read data/category different types/merge data/rename data
% Read data
data = readtable('hurricanes90s.txt');
data.Country = categorical(data.Country); % change 'character' into character;

% TODO - Add a new variable named 'Location' to the table 'data'
% It should be categorical  and have a value 'Land' if the row has a valid
% country listed and 'N/A' otherwise.

% Find the distinct categories
types = categories(data.Country);

% Remove 'N/A' from cats
countryNames = setdiff(types,'N/A');

% Merge all countries/replace countryNames with land in data.Country
data.Location = mergecats(data.Country,countryNames,'Land')

% Rename N/A to Sea
data.Location = renamecats(data.Location,'N/A','Sea')

% Find index of on-land and on-sea observation index/find the location of
% selected data (0/1)
onland = data.Location == 'Land';
onsea = data.Location == 'Sea';

avgP = mean(data.Pressure,'omitnan'); % ignore missing values

% deal with missing data
im = ismissing(data); % logic, scan all the table
mrow = any(im,2); % 1 represents by column (horizontally scan), 2 means by row (vertically scan);
data(mrow,:) = []; % delete row with missing data

% Create scatter plots/find the correspondant data (in the same row) (1-plot,0-not plot)
scatter(data.Windspeed(onsea),data.Pressure(onsea),'b','filled'); % or
scatter(data.Windspeed(~onsea),data.Pressure(~onsea),'r','filled'); % ~ means NOT
hold on
scatter(data.Windspeed(onland),data.Pressure(onland),'r','filled');
hold off

% Optional - annotate the plot
xlabel('Wind Speed')
ylabel('Pressure')
legend('Sea','Land')

% assign level to data with certain range
binedges = [0 0.2 0.6 0.9 Inf] 
bin_numbers = discretize(x,binedges) 

cnm = {'red','green','blue','black'};
y = discretize(x,binedges,'Categorical',cnm) % rename the level with correspondant character
leGreen = nnz(y <= 'green') % <=level 2, the number of 

SSscale = [0 39 74 96 111 130 157 Inf];
catnames = {'TD','TS','1','2','3','4','5'};
data.HurrCat=discretize(data.Windspeed,SSscale,'Categorical',catnames); % add a new column with windspeed level

% plot
scatter(x,y,markersize) % change the marker size for every data point
scatter(x,y,10*x,'rs','filled','MarkerFaceAlpha',0.2) % change the marker size with x/y, modify marker filling/transparency

xlim([x1,x2]);
grid('on');
grid('off');
grid('minor');
axis('tight'); 
axis('square'); 

%%% project1
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
scatter(data.CityMPG(MPGClass == 'Low'),data.HighwayMPG(MPGClass == 'Low'),'r','filled')
hold on
scatter(data.CityMPG(MPGClass == 'Medium'),data.HighwayMPG(MPGClass == 'Medium'),'b','filled')
hold on
scatter(data.CityMPG(MPGClass == 'High'),data.HighwayMPG(MPGClass == 'High'),'k','filled')
hold off
grid on
xlabel('City MPG')
ylabel('Highway MPG')
legend('Low Combined MPG','Medium Combined MPG','High Combined MPG');

% datastore
ds = datastore('dirName/fileName.txt'); % directory name (all the files) or file name (single);
dat = datastore('hurricaneData/Location/hurricaneData1990.txt');
ds; % return all the properties
ds.propertyName % use dot notation to access the property value.
v = dat.VariableNames; % the name of variables;
d = data.Delimiter; % delimiter property;
fname = dat.Files; % access the Files property (name and directory) of dat and save the results to a variable named fname

dat.CommentStyle='##'; % recognize comment, and do not import comment into data;
preview(dat) % Preview the contents of hurricaneData1990.txt using the dat datastore
ls('folder'); % list the files in the folder/directory

bds.ReadVariableNames = false; % set as false/0 if there is no heading in the data
bds.VariableNames={'name1','name2'}; % set headings manually;
bds.TextscanFormats; % format of each variables/columns 
% '%C'-categorical, '%f'-numeric, '%q'-non_numeric...
% strings,'%D'-datetime/'%T'-time 
bds.TextscanFormats{1}='%C'; % change the type of the column;

ds = datastore('airlines.csv','TreatAsMissing','NA'); % easier step;
ds=datastore('electricityData.txt','Delimiter','||'); % recognize the delimter
selected_set=ds.VariableNames; % read all the variable names;
dat.SelectedVariableNames = selected_set([1,2,4]); % the 1,2,4 columns to import from;
dat.SelectedVariableNames = {'Number','Timestamp','Country'}; % set the variables to import;
locdata=read(dat); % only import the selected variables/columns

%%% 
% 1-creat datastore,2-modify datastore properties,3-import use read function
dat = datastore('hurricaneData/Location/'); % import all files in the location into datastore;
% Adjust properties
dat.CommentStyle = '##';
% TODO - Import data
hurrs1 = read(dat);
hurrs2 = read(dat);
reset(dat);
hurrs = readall(dat)

dat.ReadSize = 15; % set lines to import data from
r1 = read(dat); % first 15 lines
r2 = read(dat); % all the rest lines

reset(dat) % but the readsize is still 15
r = readall(dat); % read all lines

% merge data
T12 = join(T1,T2);
T12 = innerjoin(T1,T4);
T12 = outerjoin(T1,T2,'MergeKeys',true);

%%% Create tables for location, windspeed, and pressure
readHurricaneData;
% TODO: Combine locdata, wsdata, and pdata into one table named hurr
hurr = join(locdata,wsdata);
hurr = outerjoin(hurr,pdata,'Mergekeys',true);

% 6.2










