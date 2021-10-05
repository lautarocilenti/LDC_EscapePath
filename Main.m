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

    %Initial Conditions

    [xoSet] = GenerateInitialConditionsFloquet(theta,M);

    %Distributed Paths
    phiSetRaw = IntegrateRHS(xoSet,M);
    % 

    phiSet = PostProcessTrajectories2(phiSetRaw,M);
    
    %Distributed Path Energies
    [S] = IntegrateLagrangian(phiSet,M);
    
    phiSet = ReduceSizeToPeriods(phiSet,M);


    %In a loop augment initial conditions near low energy local minima
    Descent.Count = 0;
    Descent.newStart = true;
    msLog = {{theta,S,Descent}};
    save("Data/ActionPlot/initialsearch.mat",'M','phiSet','S','theta','Descent','msLog','-v7.3')
else
    
     data = load('Data/ActionPlot/initialsearch.mat');
         
 
    theta = data.theta;
    S = data.S;
    
    Descent.Count = 0;
    Descent.newStart = true;
    msLog = {{theta,S,Descent}};
    phiSet = data.phiSet;
    M = data.M;
    mTemp = Parameters();
    M.MS = mTemp.MS;
    M.descent = mTemp.descent;
    
end

[phiSet,msLog] = RunMinSearch(phiSet,msLog,M); 
theta = msLog{end}{1}; S = msLog{end}{2};
% 
[minS,minPhiIndex] = IdentifyMPEP(S);

toc

%Data output
data.minPhiIndex = minPhiIndex; data.minS = minS; data.theta = theta;
data.S = S; data.phiSet = phiSet;
data.M = M; data.attractors = M.Mrhs.FixedPoints.FP;
data.msLog = msLog;

% save('temp.mat');
% SaveToFile(data,M);
if ~CheckIfCluster()
    PlotGenerator(data)
end

% CleanUpParpool();



end

