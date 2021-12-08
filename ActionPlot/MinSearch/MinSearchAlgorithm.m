function [phiSetOut,msLog,TerminateFlag] = MinSearchAlgorithm(phiSet,msLog,nLM,M)
%MINSEARCHALGORITHM             
% Runs search algorithm. Chooses initial starting points for algorith then
% it indeces whether the starting point runs a gradient method or a
% nongradient method. When running a gradient method all starting points
% run a variation of gradient descent until a discontinuity is detected.
% When a discontinuity is detected that point is logged to continue with a
% non gradient method. 


%extract data from min search log
theta = msLog{end}{1};
cost = msLog{end}{2};
searchPrev = msLog{end}{3};




if searchPrev.newStart  %first loop

    [~,cCurrent,iLM] = IdentifyInitialSearchValues(theta,cost,nLM,M);
    thetaCurrent = theta(:,iLM);
    runTheta = true(size(cCurrent));
    stepSize = M.descent.Gamma*ones(size(cCurrent));
    fdCostPrev = NaN(size(thetaCurrent));
    iRepeat = false(size(cCurrent));

    
    if contains(M.searchAlgorithm,"Gradient")
        iNonGradient = false(size(cCurrent));
        State = {};
    elseif contains(M.searchAlgorithm,"Simplex")
        State.simplexState = ones(size(cCurrent));
        State.thetaPrev = thetaCurrent;
        State.sPrev = cCurrent;
        iNonGradient = true(size(cCurrent));
    elseif contains(M.searchAlgorithm,"Fletcher-Reeves")
        iNonGradient = false(size(cCurrent));
        State.state = ones(size(cCurrent));
    else
        State = {};
        iNonGradient = true(size(cCurrent));
    end

    
else
    
    c = searchPrev.Count;
    thetaCurrent = searchPrev.thetaNew{c};
    cCurrent =  searchPrev.sNew(c,:);
    runTheta =  searchPrev.runTheta(c,:);
    fdCostPrev = searchPrev.fdCost{c};
    stepSize = searchPrev.stepSize{c};
    iRepeat = searchPrev.iRepeatNext;
    iNonGradient = searchPrev.iNonGradient;
    State = searchPrev.State;
    

    
end

[thetaNewSearch,iNewSearch,fdCost,stepDirectionSearch,State] = ... 
    DetermineNextStep(thetaCurrent,cCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,stepSize,State,M);


[cNewSearch,phiSetNewSearch] = CostFunction(thetaNewSearch,M);

thetaNew = thetaCurrent;
cNew = cCurrent;
if contains(M.searchAlgorithm,"Simplex")
    thetaNew(:,iNewSearch) = thetaNewSearch;
    cNew(iNewSearch) = cNewSearch;
else
    [iBestNewValue] = BestNewSteps(cNewSearch,iNewSearch);
    thetaNew(:,runTheta) = thetaNewSearch(:,iBestNewValue);
    cNew(runTheta) = cNewSearch(iBestNewValue);
end


[thetaOut,cOut,phiSetOut] = MergeNewData(theta,cost,phiSet,thetaNewSearch,cNewSearch,phiSetNewSearch);
   

[thetaNew,cNew,stepSize,iNonGradient,iRepeatNext,State] =  ... 
    EvaluateNewStep(thetaNew,cNew,thetaCurrent,cCurrent,iNonGradient,runTheta,stepSize,State,M);


iStop = find(stepSize<= M.descent.minGamma);
runTheta(iStop) = false;


if searchPrev.Count>=1 & sum(runTheta)>M.minICStopRemoval & ~contains(M.searchAlgorithm,"Simplex")%every iteration remove worst cost search until min
    [~,isortCostNew] = sort(cNew,'descend');
    sortedRunTheta = runTheta(isortCostNew);
    iStop2 = isortCostNew(min(find(sortedRunTheta==true)));
    runTheta(iStop2) = false;
    cNew
    runTheta
end

TerminateFlag = ~any(runTheta);
    

    if searchPrev.newStart
        searchPrev.Count = 0;
    else
        search = searchPrev;
    end
    search.newStart = false;
    search.Count = searchPrev.Count+1;
    c = search.Count;

    search.fdCost{c} = fdCost;
    search.runTheta(c,:) = runTheta;
    search.theta{c} = thetaOut;
    search.S{c} = cOut;
    search.thetaCurrent{c}= thetaCurrent;
    search.sCurrent(c,:) = cCurrent;
    search.thetaNew{c} = thetaNew;
    search.sNew (c,:) = cNew;
    search.stepSize{c} = stepSize;
    search.iRepeatNext = iRepeatNext;
    search.iNonGradient = iNonGradient;
    search.stepDirectionSearch{c} = stepDirectionSearch;
    search.iNewSearch{c} = iNewSearch;
    search.State = State;
    msLog{end} = {thetaOut,cOut,search};

end






