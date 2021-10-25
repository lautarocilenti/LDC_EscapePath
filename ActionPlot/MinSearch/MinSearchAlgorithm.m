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

    [~,sCurrent,iLM] = IdentifyInitialSearchValues(theta,S,nLM);
    thetaCurrent = theta(:,iLM);
    runTheta = true(size(sCurrent));
    stepSize = M.descent.Gamma*ones(size(thetaCurrent));
    fdCostPrev = NaN(size(thetaCurrent));
    iRepeat = false(size(sCurrent));
    
    if contains(M.searchAlgorithm,"Gradient")
        iNonGradient = false(size(sCurrent));
    else
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
    

    
end

[thetaNewSearch,iNewSearch] = ... 
    DetermineNextStep(thetaCurrent,sCurrent,runTheta,fdCostPrev,iRepeat,iNonGradient,stepSize,M)


[sNewSearch,~] = CostFunction(thetaNewSearch,M);

[iBestNewValue] = BestNewSteps(sNewSearch,iNewSearch);

thetaNew = thetaCurrent;
sNew = sCurrent;
thetaNew(:,runTheta) = thetaNewSearch(:,iBestNewValue);
sNew(runGSOnTheta) = sNewSearch(iBestNewValue);
[thetaOut,sOut] = MergeNewData(theta,s,thetaNewSearch,sNewSearch);
   

[thetaNew,sNew,stepSize,iNonGradient,iRepeatNext] =  ... 
    EvaluateNewStep(thetaNew,sNew,thetaCurrent,sCurrent,iNonGradient,runTheta,M);


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

    search.runGSOnTheta(c,:) = runGSOnTheta;
    search.theta{c} = thetaOut;
    search.S{c} = sOut;
    search.thetaCurrent{c}= thetaCurrent;
    search.sCurrent(c,:) = sCurrent;
    search.thetaNew{c} = thetaNew;
    search.sNew (c,:) = sNew;
    search.stepSize{c} = stepSize;
    search.iRepeatNext = iRepeatNext;
    search.iNonGradient = iNonGradient;
    msLog{end+1} = {thetaOut,sOut,search};

end






