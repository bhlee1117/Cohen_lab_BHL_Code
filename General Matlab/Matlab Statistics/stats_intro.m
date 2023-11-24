% Statistics Bootcamp
% JHH 120823


%% Useful hypothesis functions:
% http://www.mathworks.com/help/toolbox/stats/f31338.html

%% Statistical arrays


%% Categorical arrays: discrete data

clear all;
load fisheriris;

% Create categorical arrays
n1 = nominal(species);
n2 = nominal(species,{'species1','species2','species3'});
o1 = ordinal(n2,{},{'species1','species3','species2'});
o2 = sort(o1);

moved_up = (o1 < o2)

% Reorder levels 
labels2 = getlabels(o2)
o3 = reorderlevels(o2,labels2([1 3 2]));
labels3 = getlabels(o3)
o4 = sort(o3);

% Rename levels
o5 = o4;
o5(1:50) = 'low';
o5(51:end) = 'high';
getlabels(o5)
o5 = droplevels(o5,{'species1','species2','species3'});
getlabels(o5);

o5 = mergelevels(o4,{'species1'},'low');
o5 = mergelevels(o5,{'species2','species3'},'high');
getlabels(o5)

% Create and combine arrays
sl = meas(:,1); % Sepal length data
sw = meas(:,2); % Sepal width data
SL1 = ordinal(sl,{'short','long'},[],...
              [min(sl),median(sl),max(sl)]);
SW1 = ordinal(sw,{'short','long'},[],...
              [min(sw),median(sw),max(sw)]);
S1 = [SL1,SW1];
S1(1:10,:)          

SL2 = nominal(sl,{'short','long'},[],...
              [min(sl),median(sl),max(sl)]);
SW2 = nominal(sw,{'skinny','wide'},[],...
              [min(sw),median(sw),max(sw)]);
S2 = [SL2,SW2];
getlabels(S2)

SetosaObs = ismember(n1,'setosa');
SetosaObs = (n1 == 'setosa');
SetosaData = meas(SetosaObs,:);

%% Dataset arrays for heterogeneous data
clear all;
load fisheriris


% Create datasets
NumObs = size(meas,1);
NameObs = strcat({'Obs'},num2str((1:NumObs)','%-d'));
iris = dataset({nominal(species),'species'},...
               {meas,'SL','SW','PL','PW'},...
               'ObsNames',NameObs);
iris(1:5,:)

% Add metadata
desc = 'Fisher''s iris data (1936)';
units = [{''} repmat({'cm'},1,4)];
info = 'http://en.wikipedia.org/wiki/R.A._Fisher';
iris = set(iris,'Description',desc,...
                'Units',units,...
                'UserData',info);
get(iris)

% Access data
SepalObs = iris([1,3,5],2)
SepalObs = iris({'Obs1','Obs3','Obs5'},'SL')

BigSWLengths = iris.SL(iris.SW > 4)

varname = 'SW';
iris.(varname)(1:10)

% Join data

SepalData = iris(:,{'SL','SW'});
PetalData = iris(:,{'PL','PW'});
newiris = [SepalData,PetalData];
size(newiris)

newiris.SepalData = [newiris.SL,newiris.SW];
newiris.PetalData = [newiris.PL,newiris.PW];
newiris(:,{'SL','SW','PL','PW'}) = [];
size(newiris)

snames = nominal({'setosa';'versicolor';'virginica'});
CC = dataset({snames,'species'},{[38;108;70],'cc'})
iris2 = join(iris,CC);
iris2([1 2 51 52 101 102],:)

% Compute with datasets
summary(iris)
summary(newiris)

SepalMeans = mean(newiris.SepalData)
SepalMeans2 = datasetfun(@mean,newiris,'UniformOutput',false)

% Group functions
[order,number,group_mean,group_median,group_iqr] = grpstats(meas,species,{'gname','numel','mean',@median,@iqr})

%% Statistical visualization:

% Scatter plots: gscatter, gplotmatrix

clear all;
figure;
gscatter(iris.SL, iris.SW, iris.species, '', 'xos');
xlabel('Sepal Length')
ylabel('Sepal Width')

figure;
xvars = [iris.SL, iris.PL];
yvars = [iris.SW, iris.PW];
gplotmatrix(xvars,yvars,iris.species,'','xos')

% Box plots: boxplot
figure;
boxplot(iris.SL, iris.species, 'notch', 'on');

% Distribution plots: cdfplot, normplot, qqplot, 
clear all;
y = evrnd(0,3,100,1);
cdfplot(y)
hold on
x = -20:0.1:10;
f = evcdf(x,0,3);
figure;
plot(x,f,'m')
legend('Empirical','Theoretical','Location','NW')

figure;
x = poissrnd(10,50,1);
y = poissrnd(5,100,1);
qqplot(x,y);

figure;
x = normrnd(10,1,25,1);
normplot(x)

figure;
x1 = wblrnd(3,3,100,1);
x2 = raylrnd(3,100,1);
probplot('weibull',[x1 x2])
legend('Weibull Sample','Rayleigh Sample','Location','NW')

%% Hypothesis testing

% ztest, ttest, lillietest

clear all;
load gas
prices = [price1 price2];

% Check if from normal distribution with Lilliefours test
figure;
normplot(prices);

lillietest(price1);
lillietest(price2);

% ztest with known sigma
[h,pvalue,ci] = ztest(price1/100,1.15,0.04)

%ttest with unknown sigma
[h,pvalue,ci] = ttest(price2/100,1.15)

% ttest2 -- two indenpendent samples from normal dist with same unknown
% sigma with same mean or unequal means
[h,sig,ci] = ttest2(price1,price2)

figure;
boxplot(prices, 1)
set(gca,'XTick',[1 2])
set(gca,'XtickLabel',{'January','February'})
xlabel('Month')
ylabel('Prices ($0.01)')

%% ANOVA

% 1-way ANOVA
clear all;
load hogg
[p,tbl,stats] = anova1(hogg);
[c,m] = multcompare(stats)

% 2-way ANOVA -- must be balanced
clear all;
load popcorn;
anova2(popcorn, 3);

% n-way ANOVA
y = [52.7 57.5 45.9 44.5 53.0 57.0 45.9 44.0]';
g1 = [1 2 1 2 1 2 1 2]; 
g2 = {'hi';'hi';'lo';'lo';'hi';'hi';'lo';'lo'}; 
g3 = {'may';'may';'may';'may';'june';'june';'june';'june'};
p = anovan(y,{g1 g2 g3});

%% MANOVA

clear all;
load carsmall;
x = [MPG Horsepower Displacement Weight];

figure;
gplotmatrix(x,[],Model_Year,[],'+xo');
[d,p,stats] = manova1(x,Model_Year);
