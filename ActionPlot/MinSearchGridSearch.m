function [phiSetOut,msLog,TerminateFlag] = MinSearchGridSearch(phiSet,msLog,nLM,M)
%MINSEARCH1D                 
% Identify all local minima (all nearest neighbors have larger energy)
% Select up to nLM local minima with smallest energy
% Generate initial conditions between the neighbors to the local minima and the local minima.
% Run rise and fall method for new initial conditions


%extract data from min search log
theta = msLog{end}{1};
s = msLog{end}{2};
gsPrev = msLog{end}{3};




    if gsPrev.newStart  %first loop

        [~,sCurrent,iLM] = IdentifySmallValues(theta,s,nLM);
        thetaCurrent = theta(:,iLM);
        runGSOnTheta = true(size(sCurrent));
        stepSize = M.descent.Gamma*ones(size(thetaCurrent));
    else
        c = gsPrev.Count;
        thetaCurrent = gsPrev.thetaNew{c};
        sCurrent =  gsPrev.sNew(c,:);
        runGSOnTheta =  gsPrev.runGSOnTheta(c,:);
        stepSize = gsPrev.stepSize{c};

    end
%     PL_ThetaSpace2(theta,s)
%     PL_ThetaSpace3(thetaCurrent)
%     view(0,0)
%     drawnow()
%     hold off
%     dummy = 1;
    [thetaNewGridSearch,iGridSearch] = GridSearch(theta,thetaCurrent(:,runGSOnTheta),find(runGSOnTheta),stepSize,M);


    thetaNewGridSearch(3,:) = mod(thetaNewGridSearch(3,:),2*pi);
    thetaNewGridSearch(1:2,:) = mod(thetaNewGridSearch(1:2,:),pi);

    
    [sNewGridSearch,phiSetNewGridSearch] = CostFunction(thetaNewGridSearch,M);


    iGridSearchUnique = unique(iGridSearch);
    iNewValue = zeros(1,length(iGridSearchUnique));
    for i = 1:length(iGridSearchUnique)
        ii = find(iGridSearchUnique(i) == iGridSearch);
        [~,imin] = min(sNewGridSearch(ii));
        iNewValue(i) = ii(imin);
    end
    
    thetaNew = thetaCurrent;
    sNew = sCurrent;
    thetaNew(:,runGSOnTheta) = thetaNewGridSearch(:,iNewValue);
    sNew(runGSOnTheta) = sNewGridSearch(iNewValue);
    [thetaOut,phiSetOut,sOut] = MergeNewData(theta,phiSet,s,thetaNewGridSearch,phiSetNewGridSearch,sNewGridSearch);
   


    

    iCancelMove = (sNew>=sCurrent);
    thetaNew(:,iCancelMove) = thetaCurrent(:,iCancelMove);
    sNew(iCancelMove) = sCurrent(iCancelMove);
    stepSize(iCancelMove) = stepSize(iCancelMove)/1.1;
    
    iStop = find(stepSize<= 1E-5);
    runGSOnTheta(iStop) = false;
    
    TerminateFlag = ~any(runGSOnTheta);
    

    if gsPrev.newStart
        gsPrev.Count = 0;
    else
        gs = gsPrev;
    end
    gs.newStart = false;
    gs.Count = gsPrev.Count+1;
    c = gs.Count;

    gs.runGSOnTheta(c,:) = runGSOnTheta;
    gs.theta{c} = thetaOut;
    gs.S{c} = sOut;
    gs.thetaCurrent{c}= thetaCurrent;
    gs.sCurrent(c,:) = sCurrent;
    gs.thetaNew{c} = thetaNew;
    gs.sNew (c,:) = sNew;
    gs.stepSize{c} = stepSize;
    msLog{end+1} = {thetaOut,sOut,gs};

end

function [thetaOut,phiSetOut,sOut] = MergeNewData(theta,phiSet,S,thetaNew,phiSetNew,sNew)
    thetaOut = [theta thetaNew];
    phiSetOut = [phiSet phiSetNew];
    sOut = [S  sNew];

end

function [thetaOut,sOut,iOut] = IdentifySmallValues(theta,S,nLM)

    [~,is] = sort(S,'ascend');
    if length(is) > nLM
        is = is(1:nLM);
    end
    thetaOut = theta(:,is);
    sOut = S(is);
    iOut = is;
end

%

function [thetaNew,iGridSearch] = GridSearch(theta,thetaCurrent,iCurrent,stepSize,M)
   
    d = size(thetaCurrent,1);
    n = size(thetaCurrent,2)
    nNewPoints = d*2;
    thetaNew = zeros(d,nNewPoints*n);
    iGridSearch = zeros(1,nNewPoints*n);
   
    for i = 1:n
        j = (i-1)*(nNewPoints)+1;
        k = (i)*(nNewPoints);
        v1 = rand(1,d);
        v = [v1'/norm(v1) null(v1)];
        step = stepSize(i)*[v -v];
        thetaNew(:,j:k) = thetaCurrent(:,i)+step;
        iGridSearch(j:k) = iCurrent(i);
    end
end


    

% 
% 
% function [thetaNew,iGridSearch] = GridSearchNearDiscontinuity(theta,s,thetaNearDisc,iDiscThetaFind,sNearDisc,M)
%     d = size(thetaNearDisc,1);
%     nNewPoints = d*2-1;
%     thetaNew = zeros(d,nNewPoints*size(thetaNearDisc,2));
%     iGridSearch = zeros(1,nNewPoints*size(thetaNearDisc,2));
%     for i = 1:size(thetaNearDisc,2)
%         iNewStart = (i-1)*(nNewPoints)+1;
%         distance = vecnorm(theta-thetaNearDisc(:,i),2,1);
%         [distanceSorted,iSortDistance] = sort(distance,'ascend');
%         sSorted = s(iSortDistance);
%         thetaSorted = theta(:,iSortDistance);
%         iClosestDiscPoint = min(find(sSorted>M.descent.DiscThresh*sNearDisc(i)));
%         if isempty(iClosestDiscPoint)
%             continue
%         end
%         thetaAtNearDisc = thetaSorted(:,iClosestDiscPoint);
%         dirMagnitude = distanceSorted(iClosestDiscPoint);
%         dirVector = (thetaAtNearDisc-thetaNearDisc(:,i))./dirMagnitude;
%         orthogonalVectors = null(dirVector');
%         L = size(orthogonalVectors,2);
%         
%         thetaNew(:,iNewStart) = thetaNearDisc(:,i)+.5*dirMagnitude*dirVector;
%         thetaNew(:,iNewStart+1:iNewStart+L) = thetaNew(:,iNewStart)+.1*dirMagnitude*orthogonalVectors;
%         thetaNew(:,iNewStart+1+L:iNewStart+2*L) = thetaNew(:,iNewStart)-.1*dirMagnitude*orthogonalVectors;
%         iGridSearch(iNewStart:iNewStart+2*L) = iDiscThetaFind(i);
%     end
%         iGridSearch(iGridSearch==0) = [];
%         thetaNew(:,all(thetaNew==0,1)) = [];
%        
%         
%     
% end


function [SNew,phiSetNew] = CostFunction(theta,M)
%     M.progressbar = false;
     phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
     SNew = IntegrateLagrangian(phiSetNew,M); 
end

