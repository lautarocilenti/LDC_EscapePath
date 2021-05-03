function [minTheta,minPhi,minS] = RunMultiStart(theta,M)
%RUNGLOBALSEARCH 
eps= .5;
ms = MultiStart('UseParallel',true,'Display','iter');
f = @(theta) FindPathAndEnergy(theta,M);
opts = optimoptions(@fminunc);
problem = createOptimProblem('fminunc','x0',theta,'objective',f,'lb',0,'ub',2*pi,'options',opts);
minTheta = run(ms,problem,20);
[minPhi] =IntegrateRHS(GenerateInitialConditions(minTheta,M),M);
minS = IntegrateLagrangian(minPhi,M);
end

function [S] = FindPathAndEnergy(thetaSet,M)
    [S] = IntegrateLagrangian(IntegrateRHS(GenerateInitialConditions(thetaSet,M),M),M);
end

