function [minTheta,minPhi,minS] = RunNoConstraintOpt(x0,M)
%RUNNOCONSTRAINTOPT 
eps= .5;
f = @(theta) FindPathAndEnergy(theta,M)
options = optimoptions('fminunc','Display','iter','OptimalityTolerance',1E-12);
[minTheta,minS] = fminunc(f,x0,options)
[minPhi] =IntegrateRHS(GenerateInitialConditions(minTheta,M),M);
end

function [S] = FindPathAndEnergy(thetaSet,M)
    [S] = IntegrateLagrangian(IntegrateRHS(GenerateInitialConditions(thetaSet,M),M),M);
end
