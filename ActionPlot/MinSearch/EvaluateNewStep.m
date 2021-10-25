function [thetaNew,sNew,stepSize,iNonGradient,iRepeatNext] = EvaluateNewStep(thetaNew,sNew,thetaCurrent,sCurrent,iNongradient,runTheta,M);
%EVALUATENEWSTEP 

%nongradient cancellations
iCancelNonGradientMove = (sNew>=sCurrent & iNongradient & runTheta);


%gradient descent cancellations
iGradient = ~iNonGradient & runTheta;
iCancelGradientMove = (sNew>=M.descent.DiscThresh*sCurrent & iGradient); %indices that increased the cost by more than 20% 
iRepeatNext = false(size(sCurrent));
iRepeatNext(iCancelGradientMove) = true;
stepSize(:,iCancelGradientMove) = stepSize(:,iCancelGradientMove)/2;
iNonGradient(iCancelGradientMove) = true;

iCancelMove = iCancelGradientMove | iCancelNonGradientMove;
thetaNew(:,iCancelMove) = thetaCurrent(:,iCancelMove);
sNew(iCancelMove) = sCurrent(iCancelMove);
stepSize(iCancelMove) = stepSize(iCancelMove)/2;
    

end

