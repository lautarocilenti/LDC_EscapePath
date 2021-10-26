function [fdCost] = ApproximateGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iNonGradient,M)
%APPROXIMATEGRADIENT 
    fprintf("Approximating Gradient ... \n")
    if M.xcoordinates
        [fdCost] = FiniteDifferencePositiveGradientOnXCoordinates(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iNonGradient,M)
    else
        [fdCost] = FiniteDifferencePositiveGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iNonGradient,M);
    end
    fprintf(" ... Completed Approximating Gradient \n")
end



function [fdCost] = FiniteDifferencePositiveGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iNonGradient,M)
    if M.methodTest
        d = M.dim;
    else
        d = M.dim-1;
    end
    
    fdStep = M.descent.fdStep*ones(size(sCurrent));
    iFindFDCost = find(~iRepeat & runDescentOnTheta &~iNonGradient); %check which indices need new gradient
    fdCost = fdCostPrev; %load prior gradient values
    if any(iFindFDCost)
        fdTheta = zeros(d,d);
        fdThetaAll = zeros(d,d*length(iFindFDCost));
        fdSInit = zeros(1,d*length(iFindFDCost));
        fdStepDivisor =  zeros(1,d*length(iFindFDCost));
        for ii = 1:length(iFindFDCost)
            for id = 1:d
                iFD = iFindFDCost(ii);
                fdStepVector = zeros(d,1);
                fdStepVector(id,1) = fdStep(iFD);
                fdTheta(:,id) = thetaCurrent(:,iFD)+fdStepVector;
            end
            iA = 1+(ii-1)*d;
            iB = iA+(d-1);
            fdThetaAll(:,iA:iB) = fdTheta;
            fdSInit(1,iA:iB) = repmat([sCurrent(iFD)],1,d);
            fdStepDivisor(1,iA:iB) = repmat([fdStep(iFD)],1,d);
        end 
        [fdS,~] = CostFunction(fdThetaAll,M); %cost near initial point 
        fdGradient = (fdS - fdSInit)./fdStepDivisor;
        fdCost(:,iFindFDCost) = reshape(fdGradient,d,length(iFindFDCost));
    end
end



function [fdCost] = FiniteDifferencePositiveGradientOnXCoordinates(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iNonGradient,M)
    % x1^2+x2^2+x3^2+x4^2 = 1;
    % df1 = 2x1; df2 = 2x2; df3 = 2x3; df4 = 2x4;
    % p = point, y = coordinate
    % A = df/norm(df); %radial vector
    % B = null(A)'; orthogonal vectors to radial
    % NewX = X + fdStep*gradS
    % NewTheta = ConverToTheta(X)
    
    d = M.dim;
    dT = d - 1; %d tangent space
    
    
    fdStep = M.descent.fdStep*ones(size(sCurrent));
    iFindFDCost = find(~iRepeat & runDescentOnTheta & ~iNonGradient); %check which indices need new gradient
    fdCost = fdCostPrev; %load prior gradient values
    if any(iFindFDCost)
        fdTheta = zeros(d,d-1);
        fdThetaAll = zeros(d,dT*length(iFindFDCost));
        fdSInit = zeros(1,dT*length(iFindFDCost));
        fdStepDivisor =  zeros(1,dT*length(iFindFDCost));
        for ii = 1:length(iFindFDCost)
            iFD = iFindFDCost(ii);
            P = thetaCurrent(:,iFD);
            eRadial = (2*P)./norm(P);
            eOrthogonal = null(eRadial')';
            for id = 1:dT %iterate tangent space dimension
                fdStepVector(:,id) =  fdStep(iFD)*eOrthogonal(id,:)';
                fdTheta(:,id) = thetaCurrent(:,iFD)+fdStepVector(:,id);
            end
            iA = 1+(ii-1)*dT;
            iB = iA+(dT-1);
            fdThetaAll(:,iA:iB) = fdTheta;
            fdSInit(1,iA:iB) = repmat([sCurrent(iFD)],1,dT);
            fdOrthogonalVectors(:,iA:iB) = eOrthogonal';
            fdStepDivisor(1,iA:iB) = fdStep(iFD);
        end 
        fdThetaAllNorm = vecnorm(fdThetaAll,2,1);
        fdThetaAllNorm = repmat(fdThetaAllNorm,d,1);
        fdThetaAll = fdThetaAll./fdThetaAllNorm;
        [fdS,~] = CostFunction(fdThetaAll,M); %cost near initial point 
        fdGradientExpanded = repmat((fdS - fdSInit)./fdStepDivisor,d,1).*fdOrthogonalVectors;
        for ii = 1:length(iFindFDCost);
            iA = 1+(ii-1)*dT; iA2 = (ii)*dT;
            fdGradient(:,ii) = sum(fdGradientExpanded(:,iA:iA2),2);
        end
        fdCost(:,iFindFDCost) = fdGradient;
    end
end