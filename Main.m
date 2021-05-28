function [] = Main(parameterNames,parameterValues)
%MAIN 
close all

AddAllPaths();

%Parameters
if nargin == 2
    M = Parameters(parameterNames,parameterValues);
else
    M = Parameters();
end
tic 

%Final Condition
[M] = GetFixedPoints(M);



%Initial Conditions
theta = M.theta;
[xoSet] = GenerateInitialConditions(theta,M);

%Distributed Paths
phiSetRaw = IntegrateRHS(xoSet,M);
% 

phiSet = PostProcessTrajectories(phiSetRaw,M);
% phiSet = phiSetRaw;

%Distributed Path Energies
[S] = IntegrateLagrangian(phiSet,M);


%In a loop augment initial conditions near low energy local minima
msLog = {{theta,S}};
[phiSet,msLog] = RunMinSearch(phiSet,msLog,M); 
theta = msLog{end}{1}; S = msLog{end}{2};

[minS,minPhiIndex] = IdentifyMPEP(S);

toc

%Data output
data.minPhiIndex = minPhiIndex; data.minS = minS; data.theta = theta;
data.xoSet = xoSet; data.S = S; data.phiSet = phiSet;
data.M = M; data.attractors = M.Mrhs.FixedPoints.FP;
data.msLog = msLog;

SaveToFile(data,M);

PlotGenerator(data)





end

