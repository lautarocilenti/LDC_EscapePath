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
if strcmp(M.rhsString,"Duffing")
    M.biSet = DuffingBasins(M);
    M.Mrhs.bI = M.biSet{1}{1}; %Use first basin interpolant once a period
    M.Mrhs.A =  M.biSet{3}{1}; % use the first attractors set once a period
    M.Mrhs.r = 1E-1;
    M.Mrhs.eps = 1E-1;
else
    [basinInterpolant,attractors] = GenerateBasinInterpolant(M);
    M.Mrhs.bI = basinInterpolant;
    M.Mrhs.A = attractors;
end


%Initial Conditions
[xoSet] = GenerateInitialConditions(M.theta,M);

%Distributed Paths
phiSetRaw = IntegrateRHS(xoSet,M);
% 

[phiSet,bi,psi] = PostProcessTrajectories(phiSetRaw,M);


%Distributed Path Energies
[S] = IntegrateLagrangian(phiSet,M);

%Initial Angle for Global Search
[approxMinAngle] = FindMinAngle(S,phiSet,M);

%Search Local when eps < 2pi, global otherwise
% [minTheta,minPhi,minS] = RunGlobalSearch(approxMinAngle,M);
% [minTheta,minPhi,minS] = RunMultiStart(approxMinAngle,M);
minTheta=  approxMinAngle; minPhi = phiSet; minS = S;

toc

%Data output
data.minTheta =  minTheta; data.minPhi = minPhi; data.minS = minS;
data.xoSet = xoSet; data.S = S; data.phiSet = phiSet; 
data.basinInterpolant = M.Mrhs.bI; data.attractors = M.Mrhs.A;
data.M = M;
data.bi = bi; data.psi = psi;

SaveToFile(data,M);

PlotGenerator(data)





end

