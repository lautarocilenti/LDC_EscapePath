clear all; close all;
%first run
data = Main();

%modify initial search
theta = ConvertOToTheta(data.theta);
search.Count = 0;
search.newStart = true;
search.M = data.M;
msLog = data.msLog;
cost = data.C;
phiSet = data.phiSet;
M = data.M;
save("Data/ActionPlot/initialsearch.mat",'M','phiSet','cost','theta','search','msLog','-v7.3')

%second run
paramName = {"searchAlgorithm","continueRun","note"};
paramValue = {"Stochastic Grid",true,"SecondSearchAllTheta"};
data = Main(paramName,paramValue);




