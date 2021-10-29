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
S = msLog{end}{2};
searchPrev = msLog{end}{3};




if searchPrev.newStart  %first loop

    [~,sCurrent,iLM] = IdentifyInitialSearchValues(theta,S,nLM,M);
    thetaCurrent = theta(:,iLM);
    runTheta = true(size(sCurrent));
    stepSize = M.descent.Gamma*ones(size(sCurrent));
    fdCostPrev = NaN(size(thetaCurrent));
    iRepeat = false(size(sCurrent));

    
    if contains(M.searchAlgorithm,"Gradient")
        iNonGradient = false(size(sCurrent));
        State = {};
    elseif contains(M.searchAlgorithm,"Simplex")
        State.simplexState = ones(size(sCurrent));
        State.thetaPrev = thetaCurrent;
        State.sPrev = sCurrent;
        iNonGradient = true(size(sCurrent));
    elseif contains(M.searchAlgorithm,"Fletcher-Reeves")
        iNonGradient = false(size(sCurrent));
        State.state = ones(size(sCurrent));
    else
        State = {};
        iNonGradient = true(size(sCurrent));
    end

    
else
    
    c = searchPrev.Count;
    thetaCurrent = searchPrev.thetaNew{c};
    sCurrent =  searchPrev.sNew(c,:);
    runTheta =  searchPrev.runTheta(c,:);
    fdCostPrev = searchPrev.fdCost{c};
    stepSize = searchPrev.stepSize{c};
    iRepeat = searchPrev.iRepeatNext;
    iNonGradient = searchPrev.iNonGradient;
    State = searchPrev.State;
    

    
end

[thetaNewSearch,iNewSearch,fdCost,stepDirectionSearch,State] = ... 
    DetermineNextStep(thetaCurrent,sCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,stepSize,State,M);


[sNewSearch,phiSetNewSearch] = CostFunction(thetaNewSearch,M);

thetaNew = thetaCurrent;
sNew = sCurrent;
if contains(M.searchAlgorithm,"Simplex")
    thetaNew(:,iNewSearch) = thetaNewSearch;
    sNew(iNewSearch) = sNewSearch;
else
    [iBestNewValue] = BestNewSteps(sNewSearch,iNewSearch);
    thetaNew(:,runTheta) = thetaNewSearch(:,iBestNewValue);
    sNew(runTheta) = sNewSearch(iBestNewValue);
end


[thetaOut,sOut,phiSetOut] = MergeNewData(theta,S,phiSet,thetaNewSearch,sNewSearch,phiSetNewSearch);
   

[thetaNew,sNew,stepSize,iNonGradient,iRepeatNext,State] =  ... 
    EvaluateNewStep(thetaNew,sNew,thetaCurrent,sCurrent,iNonGradient,runTheta,stepSize,State,M);


iStop = find(stepSize<= 1E-5);
runTheta(iStop) = false;

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
    search.S{c} = sOut;
    search.thetaCurrent{c}= thetaCurrent;
    search.sCurrent(c,:) = sCurrent;
    search.thetaNew{c} = thetaNew;
    search.sNew (c,:) = sNew;
    search.stepSize{c} = stepSize;
    search.iRepeatNext = iRepeatNext;
    search.iNonGradient = iNonGradient;
    search.stepDirectionSearch{c} = stepDirectionSearch;
    search.iNewSearch{c} = iNewSearch;
    search.State = State;
    msLog{end+1} = {thetaOut,sOut,search};

end






