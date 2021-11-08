function [data] = Main(parameterNames,parameterValues)
%MAIN 
close all

AddAllPaths();
CreateParpool();

%Parameters
if nargin == 2
    M = Parameters(parameterNames,parameterValues);
else
    M = Parameters();
end
tic 

if ~M.continueRun
    %Final Condition
    [M] = GetFixedPoints(M);

    theta = GetInitialTheta(M);

    %Initial Conditions cost
    [cost,phiSet] = CostFunction(theta,M);

    


    %In a loop augment initial conditions near low energy local minima
    search.Count = 0;
    search.newStart = true;
    search.M = M;
    msLog = {};
    save("Data/ActionPlot/initialsearch.mat",'M','phiSet','cost','theta','search','msLog','-v7.3')
    msLog = {{theta,cost,search}};
else
    
     data = load('Data/ActionPlot/initialsearch.mat');
         
 
    theta = data.theta;
    cost = data.cost;
    
    search.Count = 0;
    search.newStart = true;
    msLog = data.msLog;
    
    phiSet = data.phiSet;
    M = data.M;
    mTemp = Parameters();
    M.MS = mTemp.MS;
    M.descent = mTemp.descent;
    M.searchAlgorithm = mTemp.searchAlgorithm;
    search.M = M;
    msLog{end+1} = {theta,cost,search};
    
end

[phiSet,msLog] = RunMinSearch(phiSet,msLog,M); 
theta = msLog{end}{1}; cost = msLog{end}{2};
% 
[minS,minPhiIndex] = IdentifyMPEP(cost); 

toc

%Data output
data.minPhiIndex = minPhiIndex; data.minS = minS; data.theta = theta;
data.S = ExtractAction(phiSet); data.phiSet = phiSet;
data.C = cost;
data.M = M; data.attractors = M.Mrhs.FixedPoints.FP; 
data.msLog = msLog;

% save('temp.mat');
save("Data/ActionPlot/initialsearch.mat",'M','phiSet','cost','theta','search','msLog','-v7.3')
SaveToFile(data,M);
if ~CheckIfCluster()
    PlotGenerator(data)
end

% CleanUpParpool();



end

