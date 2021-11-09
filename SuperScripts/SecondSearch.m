clear all; close all;
%first run
data = load('F:\OneDrive - University of Maryland\University of Maryland\ResearchData\ActionPlotData\November2021\CoupledDuffing\iA4\ActionPlot_TwoDuffing_Two_Oscillator_StochasticGrid_iA4_2021_11_05_14_35_32_564.mat');
data = data.data;

%modify initial search
theta = data.theta;
search.Count = 0;
search.newStart = true;
search.M = data.M;
msLog = data.msLog;
cost = data.C;
phiSet = data.phiSet;
M = data.M;
save("Data/ActionPlot/initialsearch.mat",'M','phiSet','cost','theta','search','msLog','-v7.3')

%second run
paramName = {"searchAlgorithm","continueRun","note","costType"};
paramValue = {"Stochastic Grid",true,"SecondSearchAllTheta","distancefromSaddle"};
data = Main(paramName,paramValue);




