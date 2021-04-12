function [minTheta,minPhi,minS] = RunGlobalSearch(theta,M)
%RUNGLOBALSEARCH 
eps= .5;
gs = GlobalSearch('Display','iter','NumStageOnePoints',M.gsInTrials,'NumTrialPoints',M.gsTrials);
f = @(theta) FindPathAndEnergy(theta,M);
opts = optimoptions(@fmincon,'Algorithm','sqp');
problem = createOptimProblem('fmincon','x0',theta,'objective',f,'lb',theta-eps,'ub',theta+eps,'options',opts);
minTheta = run(gs,problem);
[minPhi] =IntegrateRHS(GenerateInitialConditions(minTheta,M),M);
minS = IntegrateLagrangian(minPhi,M);
end

function [S] = FindPathAndEnergy(thetaSet,M)
    [S] = IntegrateLagrangian(IntegrateRHS(GenerateInitialConditions(thetaSet,M),M),M);
end

