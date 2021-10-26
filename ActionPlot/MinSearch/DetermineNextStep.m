function [thetaNew,iNewSearch,fdCost,stepDirectionSearch] = DetermineNextStep(thetaCurrent,sCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,stepSize,M)
%DETERMINENEXTSTEP
    iGradient = ~iNonGradient & runTheta;
    iNonGradient = iNonGradient & runTheta;
    thetaNew = [];
    iNewSearch = [];
    stepDirectionSearch = [];
    fdCost = fdCostPrev;
    if any(iGradient);
        [fdCost] = ApproximateGradient(thetaCurrent,sCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,M);
        thetaNew =  thetaCurrent(:,iGradient)-stepSize(iGradient).*fdCost(:,iGradient); %calculate new descent cost
        iNewSearch = find(iGradient);
        stepDirectionSearch = -stepSize(iGradient).*fdCost(:,iGradient);
    end
    
    if any(iNonGradient);
        [thetaNewNonGradient,iNonGradientSearch,stepDirection] = DetermineNewNonGradientPoints(thetaCurrent(:,iNonGradient),find(iNonGradient),stepSize(iNonGradient),M);
        i1 = size(thetaNew,2)+1; i2 = (i1-1) + size(thetaNewNonGradient,2);
        thetaNew(:,i1:i2) = thetaNewNonGradient;
        stepDirectionSearch(:,i1:i2) = stepDirection;
        iNewSearch(i1:i2) = iNonGradientSearch;
    end
    



end

