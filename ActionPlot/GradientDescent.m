function [theta,s,disc,descent] = GradientDescent(theta,s,m)
    d = m.dim-1;
    thetaInit = theta;
    sInit = s;
    
    %finite difference gradient
    [fdGradient] = FiniteDifferencePositiveGradient(thetaInit,sInit,M);
    
    %Stop Check
    if abs(vecnorm(fdCost,2,1))<1E-2 | abs(vecnorm(fdCost,2,1))>1E100;
        return
    end
    
    thetaNew =  thetaInit-descentStep.*fdGradient; 
    

end

function [fdGradient] = FiniteDifferencePositiveGradient(thetaInit,sInit,M)
    thetaInit = theta;
    sInit = s;
    fdStep = m.descent.fdStep;
    fdTheta = zeros(d,d);
    for id = 1:d
        fdStepVector = zeros(d,1);
        fdStepVector(id,1) = fdStep;
        fdTheta(:,id) = thetaInit+fdStepVector;
    end
    [fdS,~] = CostFunction(fdTheta,M);
    fdGradient = (fdS - sInit)./fdStep;
end