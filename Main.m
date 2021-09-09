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

    %In a loop augment initial conditions near low energy local minima
    Descent.Count = 0;
    Descent.newStart = true;
    msLog = {{theta,S,Descent}};
    save('initialSearch.mat')
else
    
    data = load("Data/Continue/continue.mat");
        

    theta = data.msLog
     msLog = data.msLog;
     phiSet = data.phiSet;
     M = data.M;
     msLog{end}{3}.newStart = true;
    M.MS.nLM = 1;
    M.MS.maxIter = 0;
    
end
% 
% load('initialSearch.mat')
%     M.MS.nLM = 3;
%     M.MS.maxIter = 1;

% Descent.Count = 0;
% Descent.newStart = true;
% msLog = {{theta,S,Descent}};
% [phiSet,msLog] = RunMinSearch(phiSet,msLog,M); 
theta = msLog{end}{1}; S = msLog{end}{2};
% 
[minS,minPhiIndex] = IdentifyMPEP(S);

toc

%Data output
data.minPhiIndex = minPhiIndex; data.minS = minS; data.theta = theta;
data.S = S; data.phiSet = phiSet;
data.M = M; data.attractors = M.Mrhs.FixedPoints.FP;
data.msLog = msLog;

% SaveToFile(data,M);

PlotGenerator(data)

CleanUpParpool();



end

