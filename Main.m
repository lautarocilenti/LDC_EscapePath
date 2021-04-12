function [] = Main(parameterNames,parameterValues)
%MAIN 
close all
tic
AddAllPaths();

%Parameters
if nargin == 2
    M = Parameters(parameterNames,parameterValues);
else
    M = Parameters();
end

%Final Condition
[basinInterpolant,attractors] = GenerateBasinInterpolant(M);
M.Mrhs.bI = basinInterpolant;

%Initial Conditions
[xoSet] = GenerateInitialConditions(M.theta,M);

%Distributed Paths
[phiSet] = IntegrateRHS(xoSet,M);

%Distributed Path Energies
[S] = IntegrateLagrangian(phiSet,M);

%Initial Angle for Global Search
[approxMinAngle] = FindMinAngle(S,phiSet,M);

%Search Local when eps < 2pi, global otherwise
[minTheta,minPhi,minS] = RunGlobalSearch(approxMinAngle,M);

toc

%Data output
data.minTheta =  minTheta; data.minPhi = minPhi; data.minS = minS;
data.xoSet = xoSet; data.S = S; data.phiSet = phiSet; 
data.basinInterpolant = basinInterpolant; data.attractors = attractors;
data.M = M;

SaveToFile(data,M);

%Plots
f1 = figure(1);
PL_Attractors(attractors);
PL_BasinBoundary(basinInterpolant);
PL_Paths(phiSet,M)

f2 = figure(2);
PL_PathEnergy(S,M)

% 
f3 = figure(3);
PL_MPEP(minTheta,minPhi,minS,M);
PL_BasinBoundary(basinInterpolant);
PL_Attractors(attractors);






end

