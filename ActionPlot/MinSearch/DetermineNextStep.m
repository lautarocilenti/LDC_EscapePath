function [thetaNew,iNewSearch,fdCost,stepDirectionSearch] = DetermineNextStep(thetaCurrent,sCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,stepSize,State,M)
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
        stepDirectionSearch = -stepSize(iGradient).*fdCost(:,iGradient);
        if M.xcoordinates
            thetaNew = thetaNew./vecnorm(thetaNew,2,1);
            stepDirectionSearch = stepDirectionSearch./vecnorm(stepDirectionSearch,2,1).*stepSize(iGradient);
        end
        iNewSearch = find(iGradient);
        
    end
    
    if any(iNonGradient);
        [thetaNewNonGradient,iNonGradientSearch,stepDirection] = DetermineNewNonGradientPoints(thetaCurrent(:,iNonGradient),sCurrent(:,iNonGradient),find(iNonGradient),stepSize(iNonGradient),State,M);
        i1 = size(thetaNew,2)+1; i2 = (i1-1) + size(thetaNewNonGradient,2);
        thetaNew(:,i1:i2) = thetaNewNonGradient;
        stepDirectionSearch(:,i1:i2) = stepDirection;
        iNewSearch(i1:i2) = iNonGradientSearch;
    end
    



end

