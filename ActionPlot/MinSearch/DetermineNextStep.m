function [thetaNew,iNewSearch] = DetermineNextStep(thetaCurrent,sCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,descentStep,M)
%DETERMINENEXTSTEP
    iGradient = ~iNonGradient & runTheta;
    iNonGradient = iNonGradient & runTheta;
    thetaNew = [];
    iNewSearch = [];
    if any(iGradient);
        [fdCost] = ApproximateGradient(thetaCurrent,sCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,descentStep,M);
        thetaNew =  thetaCurrent(:,iGradient)-descentStep(:,iGradient).*fdCost(:,iGradient); %calculate new descent cost
        iNewSearch = find(iGradient);
    end
    
    if any(iNonGradient);
        [thetaNewNonGradient,iNonGradientSearch] = DetermineNewNonGradientPoints(thetaCurrent(:,iNonGradient),find(iNonGradient),stepSize,M)
        i1 = size(thetaNew,2)+1; i2 = i1 + size(thetaNewNonGradient,2);
        thetaNew(:,i1:i2) = thetaNewNonGradient;
        iNewSearch(i1:i2) = iNonGradientSearch;
    end
    



end

