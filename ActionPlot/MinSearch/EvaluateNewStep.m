function [thetaNew,sNew,stepSize,iNonGradient,iRepeatNext,State] = EvaluateNewStep(thetaNew,sNew,thetaCurrent,sCurrent,iNonGradient,runTheta,stepSize,State,M);
%EVALUATENEWSTEP 

%EvaluateNelderMead 
if contains(M.searchAlgorithm,"Simplex")
    [thetaNew,sNew,State] = EvaluateNedlerMead(thetaNew,sNew,thetaCurrent,sCurrent,State,M);
    iRepeatNext = false(size(sCurrent));
    return
elseif contains(M.searchAlgorithm,"Fletcher-Reeves")
    iGradient = ~iNonGradient & runTheta;
    [thetaNew,sNew,State] = EvaluateFletcherReeves(thetaCurrent,thetaNew,sCurrent,sNew,iGradient,State,M);
end

%nongradient cancellations
iCancelNonGradientMove = (sNew>=sCurrent & iNonGradient & runTheta);
stepSize(iCancelNonGradientMove) = stepSize(iCancelNonGradientMove)/2;

%gradient descent cancellations
iGradient = ~iNonGradient & runTheta;
iCancelGradientMove = (sNew>=M.descent.DiscThresh*sCurrent & iGradient); %indices that increased the cost by more than 20% 
iRepeatNext = false(size(sCurrent));
iRepeatNext(iCancelGradientMove) = true;
stepSize(:,iCancelGradientMove) = M.descent.discGamma;
iNonGradient(iCancelGradientMove) = true;

iCancelMove = iCancelGradientMove | iCancelNonGradientMove;
thetaNew(:,iCancelMove) = thetaCurrent(:,iCancelMove);
sNew(iCancelMove) = sCurrent(iCancelMove);

    

end

